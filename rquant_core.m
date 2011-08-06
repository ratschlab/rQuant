function save_fname = rquant_core(CFG)
% RQUANT_CORE   Core function of rQuant.
%
%   save_fname = rquant_core(CFG)
%
%   -- input --
%   CFG:        configuration struct
%
%   -- output --
%   save_fname: file name of saved rQuant result
%
%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 3 of the License, or
%   (at your option) any later version.
%
%   Written (W) 2009-2011 Regina Bohnert, Gunnar Raetsch
%   Copyright (C) 2009-2011 Max Planck Society
%


%%%% load genes and pre-processing %%%%%
load(CFG.gene_fn, '*genes');
if exist('all_genes', 'var')
  genes = all_genes;
  clear all_genes;
end
% add eidx, adapt to closed intervals
[genes num_del num_merged] = sanitise_genes(genes, CFG);
if CFG.VERBOSE>0, fprintf(1, '\nUsing %i genes (merged %i, deleted %i)\n\n', length(genes), num_merged, num_del); end
clear num_del num_merged;
% determine ranges of transcript length bins
tl = [genes.transcript_length];
tlr = ceil([0 prctile(tl,20) prctile(tl,40) prctile(tl,60) prctile(tl,80) inf]);
CFG.transcript_len_ranges = round([tlr(1:end-1)'+1 tlr(2:end)']);
% assign length bin to each transcript
for g = 1:length(genes),
  for t = 1:length(genes(g).transcripts),
    genes(g).transcript_len_bin(t) = find(CFG.transcript_len_ranges(:,1) <= genes(g).transcript_length(t) & ...
                                          CFG.transcript_len_ranges(:,2) >= genes(g).transcript_length(t));
  end
end
if CFG.VERBOSE>1
  CFG.transcript_len_ranges
  CFG
end


%%%%% initialise profiles %%%%%
if CFG.load_profiles
  if CFG.VERBOSE>1
    fprintf(1, '\nLoading profiles... (%s)\n', CFG.profiles_fn);
  elseif CFG.VERBOSE==1
    fprintf(1, '\nLoading profiles...\n');
  end
  if strcmp(CFG.profiles_fn(end-3:end), '.mat')
    load(CFG.profiles_fn, 'profile_weights');
  else
    profile_weights = read_density_model(CFG.profiles_fn);
  end
  if CFG.norm_seqbias
    load(CFG.profiles_fn, 'seq_weights');
  end
  if ~(size(profile_weights,1)==CFG.num_plifs && size(profile_weights,2)==size(CFG.transcript_len_ranges,1))
    error('Profiles have wrong dimensions.');
  end
else
  profile_weights = ones(CFG.num_plifs, size(CFG.transcript_len_ranges,1));
  seq_weights = [];
end


%%%%% read density estimation %%%%%
if CFG.learn_profiles>0
  fprintf(1, '\nDetermining profile');
  if CFG.norm_seqbias
    fprintf(1, ' and sequence bias');
  end
  if CFG.learn_profiles==1
    fprintf(1, ' by empirical estimation...\n\n');
    assert(~CFG.norm_seqbias);
    profile_genes = genes;
    profile_weights = get_empirical_profiles(CFG, profile_genes);
  elseif CFG.learn_profiles==2
    fprintf(1, ' by optimisation...\n\n');
    if 0
    profile_genes = genes([genes.is_alt]==0); % only single-transcript genes
    tscp_len_bin = zeros(1, length(profile_genes));
    for g = 1:length(profile_genes),
      ridx = randperm(length(profile_genes(g).transcripts));
      tscp_len_bin(g) = profile_genes(g).transcript_len_bin(ridx(1));
      [tmp_coverage] = get_coverage_per_read(CFG, profile_genes(g), 1);    
      profile_genes(g).mean_ec = mean(tmp_coverage);
      if mean(tmp_coverage) < 15
        tscp_len_bin(g) = 0;
      end
    end
    % minimal number of examples in one transcript length bin
    min_bin_num = inf;
    for r = 1:size(CFG.transcript_len_ranges,1),
      min_bin_num = min(min_bin_num, sum(tscp_len_bin==r));
    end
    % take equal number of examples in each bin
    profile_genes_idx = [];
    for r = 1:size(CFG.transcript_len_ranges,1),
      fidx = find(tscp_len_bin==r);
      ridx = randperm(length(fidx));
      profile_genes_idx = [profile_genes_idx, fidx(ridx(1:min_bin_num))];
    end
    profile_genes = profile_genes(profile_genes_idx);
    %ridx = randperm(length(profile_genes));
    %num_exm = min(length(profile_genes), 500);
    %profile_genes = profile_genes(ridx(1:num_exm));
    else
      load('/fml/ag-raetsch/share/projects/rquant/data_sim/elegans/WS200/rquant/profile_genes.mat', 'profile_genes');
      %load('~/tmp/profiles.mat', 'profile_genes');
    end
    if CFG.VERBOSE>0, fprintf(1, 'Using %i genes for profile learning\n', length(profile_genes)); end
    [profile_weights, obj, seq_weights] = opt_density(CFG, profile_genes, profile_weights);
    CFG.VERBOSE = tmp_VERBOSE;
  end
    save_fname = sprintf('%s/profiles.mat', CFG.out_dir);
    if CFG.norm_seqbias
      save(save_fname, 'CFG', 'profile_genes', 'profile_weights', 'seq_weights');
    else
      save(save_fname, 'CFG', 'profile_genes', 'profile_weights');
    end
end
if CFG.VERBOSE>1
  profile_weights
end


%%%%% transcript weight estimation %%%%%
if CFG.VERBOSE>0, fprintf(1, '\nDetermining transcript weights...\n'); end
PAR.CFG = CFG;
PAR.profile_weights = profile_weights;
if CFG.norm_seqbias
  PAR.seq_weights = seq_weights;
end
if CFG.use_rproc
  num_genes_per_job = ceil(length(genes)/CFG.rproc_num_jobs);
  JOB_INFO = rproc_empty(CFG.rproc_num_jobs);
  % submit jobs to cluster
  for f = 1:CFG.rproc_num_jobs,
    PAR.genes = genes( (f-1)*num_genes_per_job+1 : min(f*num_genes_per_job, length(genes)) );
    switch CFG.organism,
     case 'elegans'
      FG.rproc_memreq = 3000;
      CFG.rproc_par.mem_req_resubmit = [6000 8000 20000];
     case 'drosophila',
      CFG.rproc_memreq = 4000;
      CFG.rproc_par.mem_req_resubmit = [7000 10000 20000];
     case 'human',
      CFG.rproc_memreq = 5000;
      CFG.rproc_par.mem_req_resubmit = [8000 12000 20000];
    end
    CFG.rproc_par.identifier = sprintf('rq.%s.f%i-', CFG.organism(1:3), f);
    fprintf(1, 'Submitting job %i (%s) to cluster\n', f, CFG.rproc_par.identifier);
    JOB_INFO(f) = rproc('opt_transcripts_caller', PAR, CFG.rproc_memreq, CFG.rproc_par, CFG.rproc_time); pause(1);
  end
  save(sprintf('~/tmp/%s_%s_%s.mat', CFG.exp, CFG.gene_source, datestr(now,'yyyy-mm-dd_HHhMM') ), 'JOB_INFO');
  % wait for jobs
  [JOB_INFO num_crashed] = rproc_wait(JOB_INFO, 60, 1, -1);
  if num_crashed>0, pause(60); end; % some delay to wait until results are written
  try
    % collect results
    for f = 1:CFG.rproc_num_jobs,
      RR = rproc_result(JOB_INFO(f), 10);
      if f>1,
        fields_genes = fieldnames(genes_new);
        fields_RR = fieldnames(RR);
        if ~isequal(fields_genes, fields_RR),
          assert(isequal(sort(fields_genes), sort(fields_RR)));
          [tmp idx1 idx2] = intersect(fields_genes, fields_RR);
          [tmp perm_idx] = sort(idx1);
          RR = orderfields(RR, idx2(perm_idx));
          fields_RR = fieldnames(RR);
          assert(isequal(fields_genes, fields_RR));
        end
      end
      genes_new((f-1)*num_genes_per_job+1 : min(f*num_genes_per_job, length(genes))) = RR;
    end
    assert(length(genes)==length(genes_new));
    genes = genes_new;
    clear genes_new;
  catch
    error('%i jobs crashed', num_crashed);
  end
else
  PAR.genes = genes;
  genes = opt_transcripts_caller(PAR);
end
% save results
save_fname = CFG.output_file;
if ~isempty(save_fname)
  if CFG.norm_seqbias
    save(save_fname, 'CFG', 'genes', 'profile_weights', 'seq_weights');
  else
    save(save_fname, 'CFG', 'genes', 'profile_weights');
  end
end

if CFG.VERBOSE>0, fprintf(1, '\n\nFinished transcript quantitation.\n'); end
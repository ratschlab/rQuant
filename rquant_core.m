%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% Written (W) 2009-2010 Regina Bohnert, Gunnar Raetsch
% Copyright (C) 2009-2010 Max Planck Society
%

function save_fname = rquant_core(CFG)
% save_fname = rquant_core(CFG)
%
% -- input --
% CFG: configuration struct
%
% -- output --
% save_fname: file name of saved rQuant result


if CFG.use_rproc,
  DEBUG = 0;
else
  DEBUG = 1;
  CFG.VERBOSE = 2;
end

%%%% load genes
load(CFG.gene_fn, '*genes');
if exist('all_genes', 'var')
  genes = all_genes;
  clear all_genes;
end
genes = genes(1:1000);

% add eidx, adapt to closed intervals
[genes num_del] = sanitise_genes(genes, CFG);

% determine ranges of transcript length bins
if isequal(CFG.gene_source, 'annotation')
  tl = [genes.transcript_length];
  if ~isequal(CFG.organism, 'human')
    tlr = ceil([0 prctile(tl,20) prctile(tl,40) prctile(tl,60) prctile(tl,80) inf])
  else
    tlr = ceil([0 prctile(tl,10) prctile(tl,20) prctile(tl,30) ...
                prctile(tl,40) prctile(tl,50) prctile(tl,60) ...
                prctile(tl,70) prctile(tl,80) prctile(tl,90) inf])
  end
  CFG.transcript_len_ranges = round([tlr(1:end-1)'+1 tlr(2:end)']);
end
% assign length bin to each transcript
for g = 1:length(genes),
  for t = 1:length(genes(g).transcripts),
    genes(g).transcript_len_bin(t) = find(CFG.transcript_len_ranges(:,1) <= genes(g).transcript_length(t) & ...
                                          CFG.transcript_len_ranges(:,2) >= genes(g).transcript_length(t));
  end
end
CFG.transcript_len_ranges


CFG
fprintf(1,'using %i genes (deleted %i)\n\n', length(genes), num_del);
clear num_del;

% initialise profiles and intron distances
if CFG.load_profiles,
  fprintf(1, '\nLoading profiles... (%s)\n', CFG.profiles_fn);
  load(CFG.profiles_fn, 'profile_weights');
  if CFG.norm_seqbias
    load(CFG.profiles_fn, 'seq_weights');
  end
  %load(CFG.profiles_fn, 'RES');
  %profile_weights = RES{end-1}.profile_weights;
  if ~(size(profile_weights,1)==CFG.num_plifs && size(profile_weights,2)==size(CFG.transcript_len_ranges,1))
    error('profiles have wrong dimensions');
  end
  x = linspace(0, 0.5, CFG.num_intron_plifs);
  intron_dists = 1-(1-x(end:-1:1))'*(1-x(end:-1:1));
  if 0
  intron_dists = RES{end-1}.intron_dists;
  if ~(size(intron_dists,1)==CFG.num_intron_plifs && size(intron_dists,2)==CFG.num_intron_plifs)
    error('intron distance matrix has wrong dimensions');
  end
  clear RES;
  end
else
  profile_weights = ones(CFG.num_plifs, size(CFG.transcript_len_ranges,1));
  x = linspace(0, 0.5, CFG.num_intron_plifs);
  intron_dists = 1-(1-x(end:-1:1))'*(1-x(end:-1:1));  
end
profile_weights(profile_weights(:,end)==0, end) = 1; % fixes 0-weights in last bin that were not learned because of too little data
profile_weights


if CFG.learn_profiles & ~CFG.load_profiles
  fprintf(1, '\nDetermining profile');
  if CFG.norm_seqbias
    fprintf(1, ' and sequence bias...\n\n');
  else
    fprintf(1, '...\n\n');
  end
  profile_genes = genes([genes.is_alt]==0); % only single-transcript genes
  tscp_len_bin = zeros(1, length(genes));
  for g = 1:length(profile_genes),
    ridx = randperm(length(profile_genes(g).transcripts));
    tscp_len_bin(g) = profile_genes(g).transcript_len_bin(ridx(1));
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
  ridx = randperm(length(profile_genes));
  num_exm = min(length(profile_genes), 500);
  profile_genes = profile_genes(ridx(1:num_exm));
  fprintf('using %i genes for profile learning\n', length(profile_genes));
  [profile_weights, obj, seq_weights] = opt_density(CFG, profile_genes);
  %[profile_weights, obj, seq_weights] = opt_density_smo(CFG, profile_genes);
  save_fname = sprintf('%s/profiles.mat', CFG.out_dir);
  if CFG.norm_seqbias
    save(save_fname, 'CFG', 'profile_genes', 'profile_weights', 'seq_weights');
  else
    save(save_fname, 'CFG', 'profile_genes', 'profile_weights');
  end
  return
end
  
fprintf(1, '\nDetermining transcript weights...\n');

PAR.CFG = CFG;
PAR.profile_weights = profile_weights;
PAR.intron_dists = intron_dists;
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
  
if ~isfield(CFG, 'out_fn')
  save_fname = sprintf('%s%s_%s_%s.mat', CFG.out_dir, CFG.exp, CFG.gene_source, CFG.method);
else
  if CFG.paired,
    save_fname = strrep(CFG.out_fn, '.mat', '_rquant_paired.mat');
  else
    save_fname = strrep(CFG.out_fn, '.mat', '_rquant.mat');
  end
end
if CFG.norm_seqbias
  save(save_fname, 'CFG', 'genes', 'profile_weights', 'seq_weights');
else
  save(save_fname, 'CFG', 'genes', 'profile_weights');
end

fprintf(1, '\n\nFinished transcript quantitation.\n');
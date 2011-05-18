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


DEBUG = 0;  

%%%% load genes
load(CFG.gene_fn, 'genes');

% merges genes, adds eidx, initialisations
[genes parent_genes num_del num_merged] = sanitise_genes(genes, CFG);

% initialise profiles and intron distances
tl = [genes.transcript_length];
if CFG.load_profiles
  if CFG.VERBOSE>1
    fprintf(1, '\nLoading profiles... (%s)\n', CFG.profiles_fn);
  else
    fprintf(1, '\nLoading profiles...\n');
  end
  profile_weights = read_density_model(CFG.profiles_fn);
  %profile_weights = load(CFG.profiles_fn);
  CFG.num_plifs = size(profile_weights,1);
  intron_dists = read_density_model(CFG.intron_dists_fn);
  %intron_dists = load(CFG.intron_dists_fn);
  CFG.num_intron_plifs = size(intron_dists,1);
  if size(intron_dists,2)~=CFG.num_intron_plifs
    error('intron distance matrix has wrong dimensions (must be the same)');
  end
  step_size = 100/size(profile_weights,2);
  tl_x = round(step_size:step_size:step_size*(size(profile_weights,2)-1));
else
  tl_x = 10:10:90;
end
% determine ranges of transcript length bins
tlr = ceil([0 prctile(tl,tl_x) inf]);
CFG.transcript_len_ranges = round([tlr(1:end-1)'+1 tlr(2:end)']);
if ~CFG.load_profiles
  fprintf(1, '\nUsing rQuant without profiles\n');
  profile_weights = ones(CFG.num_plifs, size(CFG.transcript_len_ranges,1));
  x = linspace(0, 0.5, CFG.num_intron_plifs);
  intron_dists = 1-(1-x(end:-1:1))'*(1-x(end:-1:1));  
end
% assign length bin to each transcript
for g = 1:length(genes),
  for t = 1:length(genes(g).transcripts),
    genes(g).transcript_len_bin(t) = find(CFG.transcript_len_ranges(:,1) <= genes(g).transcript_length(t) & ...
                                          CFG.transcript_len_ranges(:,2) >= genes(g).transcript_length(t));
  end
end

if CFG.subsample,
  CFG.subsample_frac_global = str2double(sprintf('%.2f', min(1,CFG.max_num_train_exm/sum([genes.exonic_len])))); 
end

% subsample here
if CFG.subsample && isequal(CFG.gene_source, 'annotation') 
  ridx = randperm(length(genes));
  CFG.subsample_idx = ridx(1:ceil(CFG.subsample_frac_global*length(genes)));
else
  CFG.subsample_idx = 1:length(genes);
end

fprintf(1,'using %i genes (merged %i, deleted %i)\n\n', length(genes), num_merged, num_del);
clear num_del;

% iteratively optimise transcript and profile weights 
iter = 0;
while (1)
  iter = iter + 1;
  fprintf(1, '\n*** Iteration %i ***\n', iter);
  fprintf(1, '\nDetermining transcript weights...\n');
  
  PAR.CFG = CFG;
  PAR.profile_weights = profile_weights;
  PAR.intron_dists = intron_dists;
  PAR.genes = genes;
  genes = opt_transcripts_caller(PAR);

  tmp = [genes.loss]; tmp = [tmp(:)];
  loss.all = [tmp.all]; loss.exons = [tmp.exons]; loss.introns = [tmp.introns];
  
  RES{iter}.profile_weights = profile_weights;
  RES{iter}.intron_dists = intron_dists;
  
  if DEBUG, keyboard; end
  
  if isequal(CFG.gene_source, 'annotation')
    save_fname = sprintf('%s%s_%s_iter%i.mat', CFG.out_dir, CFG.gene_source, CFG.method, iter);
  else
    if CFG.paired,
      save_fname = strrep(CFG.gene_fn, '.mat', '_rquant_profile_paired.mat');
    else
      save_fname = strrep(CFG.gene_fn, '.mat', '_rquant_profile.mat');
    end
  end
  genes_tmp = genes;
  [parent_genes.transcript_weights] = deal([]);
  for g = 1:length(genes),
    gidx = unique(genes(g).parent_gene);
    for n = 1:length(gidx),
      tidx = find([genes(g).parent_gene]==gidx(n));
      assert(isequal(parent_genes(gidx(n)).transcripts, {genes(g).transcripts{tidx}}));
      parent_genes(gidx(n)).transcript_weights = genes(g).transcript_weights(tidx);
    end
  end
  genes = parent_genes;
  save(save_fname, 'CFG', 'RES', 'genes');
  genes = genes_tmp;
  clear genes_tmp;

  if iter >= CFG.max_iter
    fprintf(1, '\n\nFinished transcript quantitation.\n');
    break;
  end
  
  fprintf(1, '\nDetermining profile...\n\n');
  profile_genes = genes(CFG.subsample_idx);
  % to ensure that profile_genes.transcript_weights are not empty
  del_idx = false(1,length(profile_genes));
  for g = 1:length(profile_genes),
    if isempty(profile_genes(g).transcript_weights) || any(isnan(profile_genes(g).transcript_weights))
      del_idx(g) = true;
    end
  end
  profile_genes(del_idx) = [];
  fprintf('using %i genes for profile learning\n', length(profile_genes));
  [profile_weights, intron_dists, profile_genes] = opt_profiles(CFG, profile_genes);
  
  if iter>1,
    RES{iter-1}.genes = [];
    RES{iter-1}.profile_genes = [];
  end
  RES{iter}.profile_weights = profile_weights;
  RES{iter}.intron_dists = intron_dists;
  RES{iter}.profile_genes = profile_genes;
end

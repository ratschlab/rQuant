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
end
% assign length bin to each transcript
for g = 1:length(genes),
  for t = 1:length(genes(g).transcripts),
    genes(g).transcript_len_bin(t) = find(CFG.transcript_len_ranges(:,1) <= genes(g).transcript_length(t) & ...
                                          CFG.transcript_len_ranges(:,2) >= genes(g).transcript_length(t));
  end
end

fprintf(1,'using %i genes (merged %i, deleted %i)\n\n', length(genes), num_merged, num_del);
clear num_del;

if CFG.learn_profiles & ~CFG.load_profiles
  fprintf(1, '\nDetermining profile');
  fprintf(1, '...\n\n');
  profile_genes = genes([genes.is_alt]==0); % only single-transcript genes
  tscp_len_bin = zeros(1, length(profile_genes));
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
end
  
fprintf(1, '\nDetermining transcript weights...\n');

PAR.CFG = CFG;
PAR.profile_weights = profile_weights;
PAR.genes = genes;
genes = opt_transcripts_caller(PAR);

save_fname = sprintf('%sgenes_rquant.mat', CFG.out_dir);
save(save_fname, 'CFG', 'genes', 'profile_weights');

fprintf(1, '\n\nFinished transcript quantitation.\n');

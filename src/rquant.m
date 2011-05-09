function rquant(anno_dir, track, output_file, output_dir, load_profiles, profiles_fn, intron_dists_fn, learn_profiles, num_iter, profiles_fn_out, intron_dists_fn_out)
% RQUANT   Determines the abundance of multiple transcripts per gene locus from RNA-Seq measurements.
%
%   rquant(anno_dir, track, output_file, output_dir, load_profiles, profiles_fn, 
%          intron_dists_fn, learn_profiles, num_iter, profiles_fn_out, intron_dists_fn_out)
%
%   -- input --
%   anno_dir:        directory of genes
%   track:           name of BAM file
%   output_file:     result gff3 file 
%   output_dir:      output directory
%   load_profiles:   flag that enables loading of profiles
%   profiles_fn:     name of input file that stores profiles
%   intron_dists_fn: name of input file that stores intron distances
%   learn_profiles:  flag that enables learning of profiles
%   num_iter:        number of iterations
%   profiles_fn_out: name of output file that stores profiles 
%   intron_dists_fn: name of output file that stores intron distances
%
%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 3 of the License, or
%   (at your option) any later version.
%
%   Written (W) 2009-2010 Regina Bohnert, Gunnar Raetsch
%   Copyright (C) 2009-2010 Max Planck Society
%


if exist('load_profiles', 'var')
  if ~isnumeric(load_profiles)
    load_profiles = str2num(load_profiles);
  end
end
if exist('learn_profiles', 'var')
  if ~isnumeric(learn_profiles)
    learn_profiles = str2num(learn_profiles);
  end
end
if exist('num_iter', 'var')
  if ~isnumeric(num_iter)
    num_iter = str2num(num_iter);
  end
end

% rQuant paths
global RQUANT_PATH RQUANT_SRC_PATH RQUANT_LD_LIBRARY_PATH

% interpreter paths
global INTERPRETER MATLAB_BIN_PATH OCTAVE_BIN_PATH

% optimizer related paths
global OPTIMIZER OPTIMIZER_PATH OPTIMIZER_TOOLBOX_PATH

% SAMTools path
global SAMTOOLS_DIR

addpath(sprintf('%s/mex', RQUANT_PATH));
addpath(sprintf('%s/tools', RQUANT_PATH));
addpath(sprintf('%s', RQUANT_SRC_PATH));

% more output to stdout
CFG.VERBOSE = 1;

CFG.paths = '';
CFG.out_dir = sprintf('%s/', output_dir);
unix(sprintf('mkdir -p %s', CFG.out_dir));

CFG.gene_source = 'annotation';

%%%%% directories from which to load read data and genes %%%%%
CFG.repeats_fn = '';
CFG.samtools_dir = sprintf('%s/', SAMTOOLS_DIR);

CFG.gene_dir = anno_dir; 
CFG.gene_fn = sprintf('%s/genes.mat', CFG.gene_dir);

%%%%% genome information %%%%%
try
  [s bam_header] = unix(sprintf('%s./samtools view %s -H', CFG.samtools_dir, track));
  assert(s==0);
  fidx1 = strfind(bam_header,'SN'); fidx2 = strfind(bam_header,'LN');
  for c = 1:length(fidx1), 
    contig_names{c} = bam_header(fidx1(c)+3:fidx2(c)-2); 
  end
catch
  fprintf(1, '\ncontig names could not be parsed from BAM file\n');
  return;
end

%%%%% read length, number of mapped reads, bai file %%%%%
CFG.read_len = 0;
for c = 1:length(contig_names),
  fprintf('Checking bam file: contig %i/%i\n', c, length(contig_names));
  CFG.tracks_fn{c} = {track};
  if ~exist(CFG.tracks_fn{c}{1}, 'file'),
    basedir = fileparts(CFG.tracks_fn{c}{1});
    filename = strrep(CFG.tracks_fn{c}{1}, '_files/alignments.bam', '.dat');
    unix(sprintf('mkdir %s; ln -s %s %s',  basedir, filename, CFG.tracks_fn{c}{1}));
  end 
  for f = 1:length(CFG.tracks_fn{c}),
    fname = CFG.tracks_fn{c}{f};
    if ~exist(sprintf('%s.bai', fname), 'file')
      command = sprintf('%s./samtools index %s', CFG.samtools_dir, fname);
      [s m] = unix(command);
      if ~exist(sprintf('%s.bai', fname), 'file')
        if CFG.VERBOSE>0, fprintf(1, '\nbai file for %s could not be created\n', fname); end
      end
    end
  end
  try
    [read_len mapped_reads(c)] = get_bam_properties(CFG.tracks_fn{c}{1}, CFG.samtools_dir, contig_names{c});
  catch
    read_len = 0;
    mapped_reads(c) = 0;
  end
  CFG.read_len = max(CFG.read_len, read_len);
end

%%%%% rquant parameters %%%%%
% enables taking data for both strands together
CFG.both_strands = 1;

% number of iterations (1: no profile learning)
if learn_profiles
  CFG.max_iter = num_iter;
else
  CFG.max_iter = 1;
end

% optimizer
CFG.optimizer = OPTIMIZER;

%%%%% transcript weight optimisation
% method to determine transcript weights 
CFG.method = 'pos'; % 'pos' or 'seg'
CFG.paired = 0;
% regularisation strength in transcript weight optimisation
CFG.C1 = 1;
CFG.C1_set = [0.001];
CFG.C1_loss_frac_target = 0.3;

%%%%% profile learning
% enables loading of profiles from CFG.profiles_fn
CFG.load_profiles = load_profiles;
CFG.profiles_fn = profiles_fn;
CFG.intron_dists_fn = intron_dists_fn;
% number of plifs for profile functions
CFG.num_plifs = 50;
% maximal number of positions to be considered at both transcript ends
CFG.max_side_len = 5000;
% bins for different expression levels
exr = [-1 inf];
CFG.expr_ranges = round([exr(1:end-1)'+1 exr(2:end)']);
% bins for different transcript lengths % prctile(tlen,10) prctile(tlen,90)
tlr = [0 1058 1645 2337 3382 inf];
%tlr = [0 649 1008 1379 1977 inf];
%tlr = [0 297 486 586 745 959 1451 2013 2738 4036 inf];
CFG.transcript_len_ranges = round([tlr(1:end-1)'+1 tlr(2:end)']);
% bins for distances to closest intron
CFG.num_intron_plifs = 5;
% enables subsampling of data for learning profiles
CFG.subsample = 1;
% maximal number of examples for learning profiles
CFG.max_num_train_exm = 4e6 * 5;
% fraction of genes to be subsampled for learning profiles
CFG.subsample_frac_global = 1;
% fraction of profile_genes to be subsampled for learning profiles
CFG.subsample_frac = 0.50;
% regularisation strength in profile optimisation
CFG.C2.tau   = 100;
CFG.C2.kappa = 1;
CFG.C2.theta = 10;
%CFG.C2.tau   = CFG.C2.tau;
%CFG.C2.kappa = CFG.C2.kappa;
%CFG.C2.theta = CFG.C2.theta;
CFG.C2.tau   = CFG.C2.tau*CFG.max_num_train_exm;
CFG.C2.kappa = CFG.C2.kappa*CFG.max_num_train_exm;
CFG.C2.theta = CFG.C2.theta*CFG.max_num_train_exm;

switch CFG.optimizer
 case 'mosek'
  switch determine_engine
   case 'octave'
    p = sprintf('%s/mosek/octave/:', RQUANT_PATH);
   case 'matlab'
    p = sprintf('%s/mosek/matlab/:', RQUANT_PATH);
   otherwise
    error('unknown engine %s', determine_engine);
  end
  CFG.paths = sprintf('%s:%s', p, CFG.paths);
  addpath(p);
 otherwise
  error('unknown optimizer %s', CFG.optimizer);
end

mosek_fname = which('quadprog');
if ~isempty(strfind(mosek_fname, 'toolbox/optim/optim'))
  error('mosek toolbox missing');
end


%%%%% rquant %%%%%
save_fname = rquant_core(CFG);


%%%%% write to gff file %%%%%
load(save_fname, 'genes');
fidx = find(save_fname=='/', 1, 'last');
if isempty(fidx), fidx = 0; end
unix(sprintf('ln -sf %s %s/genes_rquant.mat', save_fname(fidx+1:end), CFG.out_dir));
write_rquant_gff3(CFG, genes, output_file, mapped_reads);


%%%%% write learned read density model %%%%%
if learn_profiles,
  load(save_fname, 'RES');
  assert(length(RES)>=1);
  if CFG.max_iter > 1
    idx = length(RES)-1;
  else
    idx = 1;
  end
  write_density_model(RES{idx}.profile_weights, RES{idx}.intron_dists, profiles_fn_out, intron_dists_fn_out);
end

function rquant(anno_dir, track, output_file, output_dir, load_profiles, profiles_fn, learn_profiles, profiles_fn_out)
% RQUANT   Determines the abundance of multiple transcripts per gene locus from RNA-Seq measurements.
%
%   rquant(anno_dir, track, output_file, output_dir, load_profiles, profiles_fn, learn_profiles, profiles_fn_out)
%
%   -- input --
%   anno_dir:        directory of genes
%   track:           name of BAM file
%   output_file:     result gff3 file 
%   output_dir:      output directory
%   load_profiles:   flag that enables loading of profiles
%   profiles_fn:     name of input file that stores profiles
%   learn_profiles:  flag that enables learning of profiles
%   profiles_fn_out: name of output file that stores profiles 
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

% rQuant version
global RQUANT_VERSION

% interpreter paths
global INTERPRETER MATLAB_BIN_PATH OCTAVE_BIN_PATH

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
        if CFG.VERBOSE>0, fprintf(1, '\nbai file could not be created\n'); end
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

%%%%% transcript weight optimisation
% method to determine transcript weights 
CFG.method = 'pos'; % 'pos' or 'seg'
CFG.paired = 0;

%%%%% profile learning
% enables loading of profiles from CFG.profiles_fn
CFG.load_profiles = load_profiles;
CFG.profiles_fn = profiles_fn;
CFG.learn_profiles = learn_profiles;
% number of iterations
CFG.max_iter = 50;
% number of plifs for profile functions
CFG.num_plifs = 50;
% maximal number of positions to be considered at both transcript ends
CFG.max_side_len = 5000;
% bins for different transcript lengths % prctile(tlen,10) prctile(tlen,90)
tlr = [0 649 1008 1379 1977 inf];
%tlr = [1 750 921 1092 1263 1434 1605 1776 1947 2118 2289 2461 2632 2803 2974 3145 3316 3487 3658 3829 4000 Inf];
CFG.transcript_len_ranges = round([tlr(1:end-1)'+1 tlr(2:end)']);
% bins for distances to closest intron
CFG.num_intron_plifs = 5;
% enables subsampling of data for learning profiles
CFG.subsample = 0;
% maximal number of examples for learning profiles
CFG.max_num_train_exm = 1e6;
% fraction of positions to be subsampled for learning profiles
CFG.subsample_frac = 0.23;
% regularisation strengths
CFG.C_I = 100;
CFG.C_F = 100;
CFG.C_N = 10;


%%%%% rquant %%%%%
save_fname = rquant_core(CFG);


%%%%% write to gff file %%%%%
load(save_fname, 'genes');
fidx = find(save_fname=='/', 1, 'last');
if isempty(fidx), fidx = 0; end
write_rquant_gff3(CFG, genes, sprintf('rQuant v%s', RQUANT_VERSION), output_file, mapped_reads);


%%%%% write learned read density model %%%%%
if learn_profiles,
  load(save_fname, 'profile_weights');
  assert(length(RES)>=1);
  if CFG.max_iter > 1
    idx = length(RES)-1;
  else
    idx = 1;
  end
  write_density_model(profile_weights, profiles_fn_out);
end

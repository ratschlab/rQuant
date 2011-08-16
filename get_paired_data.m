function [pair_mat_exp, pair_mat_obs, segments, paired_reads_exp] = get_paired_data(CFG, gene, paired_reads_obs)
% GET_PAIRED_DATA   Prepares paired data for gene.
%
%   [pair_mat_exp, pair_mat_obs, segments, paired_reads_exp] = get_paired_data(CFG, gene, paired_reads_obs)
%
%   -- input --
%   CFG:              configuration struct with
%   gene:             struct defining a gene with start, stops, exons etc.
%   paired_reads_obs: struct including starts, stops and mates for observed paired-end reads
%
%   -- output --
%   pair_mat_exp:     matrix (#segments x #segments x T) of connectivity by expected paired-end reads
%   pair_mat_obs:     matrix (#segments x #segments) of connectivity by observed paired-end reads
%   segments:         distinguishable segments
%   paired_reads_exp: struct including starts, stops and mates for expected paired-end reads
%
%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 3 of the License, or
%   (at your option) any later version.
%
%   Written (W) 2011 Regina Bohnert
%   Copyright (C) 2011 Max Planck Society
%


%%% get distinguishable segments from splicegraph
segments = define_segments(gene.splicegraph{1}, gene.splicegraph{2});

%%% generate paired-end reads from annotated transcripts
ins_size = 200;
paired_reads_exp.starts = [];
paired_reads_exp.stops = [];
paired_reads_exp.mates = [];
pair_mat_exp = zeros(size(segments,1), size(segments,1), length(gene.transcripts));
offset = 0;
for t = 1:length(gene.transcripts),
  tidx = [];
  for e = 1:size(gene.exons{t},1),
    tidx = [tidx, gene.exons{t}(e,1):gene.exons{t}(e,2)];
    tidx = unique(tidx);
  end
  if length(tidx) < 2*CFG.read_len+ins_size
    continue;
  end
  cand1_starts = tidx(1:end-2*CFG.read_len-ins_size+1);
  cand1_stops = tidx(CFG.read_len-1+[1:end-2*CFG.read_len-ins_size+1]);
  cand2_starts = tidx(CFG.read_len+ins_size+[1:end-2*CFG.read_len-ins_size+1]);
  cand2_stops = tidx(2*CFG.read_len+ins_size-1+[1:end-2*CFG.read_len-ins_size+1]);
  paired_reads_exp.starts = [cand1_starts, cand2_starts];
  paired_reads_exp.stops = [cand1_stops, cand2_stops];
  paired_reads_exp.mates = [1:length(cand1_starts); length(cand1_starts)+1:length(cand1_starts)+length(cand2_starts)];
  % connectivity matrix for expected pairs
  pair_mat_exp(:,:,t) = gen_paired_segments(segments, paired_reads_exp);
  %paired_reads_exp.starts = [paired_reads_exp.starts, cand1_starts, cand2_starts];
  %paired_reads_exp.stops = [paired_reads_exp.stops, cand1_stops, cand2_stops];
  %paired_reads_exp.mates = [paired_reads_exp.mates, offset+[1:length(cand1_starts); length(cand1_starts)+1:length(cand1_starts)+length(cand2_starts)]];
  %offset = offset + length(paired_reads_exp.starts);
end

%%% connectivity  matrix for observed pairs
pair_mat_obs = gen_paired_segments(segments, paired_reads_obs);
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% Written (W) 2009-2010 Regina Bohnert, Gunnar Raetsch
% Copyright (C) 2009-2010 Max Planck Society
%

function genes = opt_transcripts_caller(PAR)
% genes = opt_transcripts_caller(PAR)
%
% -- input --
% PAR contains
%     CFG: configuration struct
%     genes: struct defining genes with start, stops, exons etc.
%     profile_weights: weights of profile functions
%     intron_dists: distances to closest intron
%
% -- output --
% genes: struct with additional fields of eg. estimated transcript weights

CFG = PAR.CFG;
genes = PAR.genes;
profile_weights = PAR.profile_weights;
intron_dists = PAR.intron_dists;
clear PAR;

%%%% paths
addpath(CFG.paths);

if CFG.VERBOSE>0, tic; end
  
switch CFG.optimizer
 case 'cplex'
  lpenv = cplex_license(0,1);
  if lpenv==0
    fprintf(1, '\n no cplex licence, using mosek\n');
    CFG.optimizer = 'mosek';
  else
    CFG.optimizer = 'cplex';
  end
 otherwise
  lpenv = 0;
end

genes(1).mean_ec = [];
genes(1).transcript_weights = [];
genes(1).transcript_weights_all = [];
genes(1).loss = struct;

chr_num = unique([genes.chr_num]);
for c = chr_num,
  chr_idx = find([genes.chr_num]==c);
  if CFG.VERBOSE>0, fprintf(1, '\nprocessing %i genes for contig %i (%s)\n', length(chr_idx), genes(chr_idx(1)).chr_num, genes(chr_idx(1)).chr); end
  for g = chr_idx,
    gene = genes(g);
    fprintf(1, '\ngene %i: %i isoform(s) with %i exonic positions\n', g, length(gene.transcripts), gene.exonic_len);
    %%%%% load exon coverage and introns for gene %%%%%
    try
      if CFG.VERBOSE>1, fprintf(1, 'Loading reads...\n'); tic; end
      if ~CFG.norm_seqbias
        [coverage excluded_reads reads_ok introns] = get_coverage_per_read(CFG, gene);
      else
        [coverage excluded_reads reads_ok introns genes(g).num_read_starts] = get_coverage_per_read(CFG, gene);
      end
      % plus strand
      fidx = find(introns(:,4)==0);
      intron_starts{1} = introns(fidx,1); 
      intron_stops{1} = introns(fidx,2);
      conf{1} = introns(fidx,3);
      % minus strand
      fidx = find(introns(:,4)==1);
      intron_starts{2} = introns(fidx,1); 
      intron_stops{2} = introns(fidx,2);
      conf{2} = introns(fidx,3);
      if CFG.VERBOSE>1, fprintf(1, 'Took %.1fs.\n', toc); end
    catch
      reads_ok = 0;
    end
    if ~reads_ok,
      if CFG.VERBOSE>0, fprintf(1, 'coverage could not be loaded for gene %i\n', g); end
      genes(g).transcript_weights(1:length(genes(g).transcripts)) = nan;
      genes(g).transcript_weights_all(1:length(genes(g).transcripts)) = nan;
      genes(g).loss.all = nan;
      genes(g).loss.exons = nan;
      genes(g).loss.introns = nan;
      continue;
    end
    genes(g).mean_ec = full(mean(sum(coverage,2)));
    genes(g).all_coverage = sum(coverage,2);
    genes(g).all_introns = introns;
    if isempty(coverage) || genes(g).mean_ec==0
      if CFG.VERBOSE>0, fprintf(1, 'no coverage for gene %i\n', g); end
      genes(g).transcript_weights(1:length(genes(g).transcripts)) = 0;
      genes(g).transcript_weights_all(1:length(genes(g).transcripts)) = 0;
      genes(g).loss.all = nan;
      genes(g).loss.exons = nan;
      genes(g).loss.introns = nan;
      continue;
    end
    %%%%% prepare intron list %%%%%
    intron_list = zeros(0,4); 
    num_transcripts(1) = length(gene.transcripts);
    if CFG.both_strands && isfield(gene, 'strands') && ~isempty(gene.strands)
      strand_str = unique(gene.strands);
      num_transcripts(2) = sum(gene.strands=='+');
      num_transcripts(3) = sum(gene.strands=='-');
    else
      strand_str = gene.strand;
      if gene.strand=='+',
        num_transcripts(2) = length(gene.transcripts);
        num_transcripts(3) = 0;
      else
        num_transcripts(2) = 0;
        num_transcripts(3) = length(gene.transcripts);
      end
    end
    assert(length(strand_str)==1|length(strand_str)==2);
    for s = 1:length(strand_str),
      strand = (strand_str(s)=='-') + 1;
      if ~isempty(intron_starts{strand})
        idx = find(intron_starts{strand}>=gene.start & intron_stops{strand}<=gene.stop);
        intron_list = [intron_list; ...
                       double([intron_starts{strand}(idx), intron_stops{strand}(idx), ...
                            conf{strand}(idx), ((strand_str(s)=='-')+1)*ones(length(idx),1)])];
      end
    end
    if CFG.VERBOSE>=2,
      fprintf(1, 'found %i transcripts, %i on + strand, %i on - strand\n', num_transcripts(1), num_transcripts(2), num_transcripts(3));
      if exist('read_starts', 'var')
        fprintf(1, 'found %i reads\n', sum(read_starts));
      end
      fprintf(1, 'found %i introns, %i on + strand, %i on - strand (%i ignored)\n', size(intron_list,1), sum(intron_list(:,4)==1), sum(intron_list(:,4)==2), size(introns,1)-size(intron_list,1));
    end
    %%%%% prepare exon and intron masks %%%%%
    exon_mask = zeros(gene.exonic_len, length(gene.transcripts));
    intron_mask = zeros(size(intron_list,1), length(gene.transcripts));
    gene.intron_dists = zeros(gene.exonic_len, length(gene.transcripts));
    for t = 1:length(gene.transcripts),
      if isfield(gene, 'strands') && ~isempty(gene.strands),
        strand_str = gene.strands(t);
      else
        strand_str = gene.strand;
      end
      dist = gen_intron_features(gene, t, CFG.num_intron_plifs, CFG.read_len);
      assert(gene.exonic_len==size(dist,1));
      idx = find(sum(sum(dist,2),3));
      for p = idx',
        tmp_dist = squeeze(dist(p,:,:));
        if all(tmp_dist(:)==0),
          gene.intron_dists(p,t) = 0;
        else
          if strand=='+'
            gene.intron_dists(p,t) = intron_dists(tmp_dist==1);
          else
            gene.intron_dists(p,t) = intron_dists(tmp_dist'==1);
          end
        end
      end
      % fill exon mask
      if strand_str=='+'
        rev_idx = 1:size(profile_weights,1);
      else
        rev_idx = size(profile_weights,1):-1:1;
      end
      exon_mask(:,t) = gen_exon_features(gene, t, CFG.num_plifs, CFG.max_side_len) * (profile_weights(rev_idx, gene.transcript_len_bin(t), gene.expr_bin(t))+1e-8) - gene.intron_dists(:,t);

      % normalise profile for sequence biases (depending on transcript sequence)
      if CFG.norm_seqbias && ~isempty(CFG.RR.seq_norm_weights),
        tidx = [];
        for e = 1:size(gene.exons{t},1),
          tidx = [tidx, gene.exons{t}(e,1):gene.exons{t}(e,2)];
          tidx = unique(tidx);
        end
        [tmp1 idx_exon_t tmp2]= intersect(gene.eidx, tidx);
        assert(isequal(tmp2,[1:length(tidx)]));
        exon_mask(idx_exon_t, t) = norm_sequence(CFG, gene, t, exon_mask(idx_exon_t, t));
      end

      % fill intron mask
      exons = gene.exons{t};
      num_found = 0;
      for e = 1:size(exons,1)-1,
        if ~isempty(intron_list)
          idx = find(exons(e,2)+1==intron_list(:,1) & exons(e+1,1)-1==intron_list(:,2));
          if CFG.both_strands && isfield(gene, 'strands') && ~isempty(idx)
            assert(all(intron_list(idx,4)==(gene.strands(t)=='-')+1));
          end
        end
        if ~isempty(intron_mask)
          if ~isempty(idx),
            num_found = num_found + 1;
          end
          intron_mask(idx,t) = 1;
        end
      end
      if num_found==0,
        if CFG.VERBOSE>0
          fprintf(1, 'introns not found in gene %i, transcript %i \n', g, t);
        end
      else
        if CFG.VERBOSE>=2
          fprintf(1, 'found %i matching introns, gene %i, transcript %i\n', num_found, g, t)
        end
      end
    end
    repeat_mask = false(gene.exonic_len, 1); 
    % fill repeat mask
    fname = sprintf('%s%s_repeat', CFG.repeats_fn, gene.chr);
    if exist(sprintf('%s.pos', fname), 'file')
      [map.pos map.repeats] = interval_query(fname, {'repeats'}, [gene.start;gene.stop]);
      if ~isempty(map.pos)
        [tmp idx1 idx2] = intersect(map.pos, gene.eidx);
        assert(length(idx2)<=length(map.pos));
        repeat_mask(idx2) = true;
      end
    end
    
    if strcmp(CFG.method, 'seg') % segment-wise
      % code only works for CFG.paired=0
      assert(size(coverage,1)==gene.exonic_len);
      segments = gen_segments(gene);
      assert(sum(segments(:,2)-segments(:,1)+1)==gene.exonic_len);
      exon_mask_pos = exon_mask;
      coverage_pos = coverage;
      exon_mask = zeros(size(segments,1),length(gene.transcripts));
      coverage = zeros(size(segments,1),1);
      for s = 1:size(segments,1),
        sidx = segments(s,1):segments(s,2);
        exon_mask(s,:) = median(exon_mask_pos(sidx,:),1);
        coverage(s) = round(median(coverage_pos(sidx)));
      end
      clear exon_mask_pos coverage_pos;
    elseif strcmp(CFG.method, 'pos') % position-wise
      % only consider regions in at least one exon
      if size(exon_mask,1)<2*(CFG.max_side_len-1)
        assert(all(any(exon_mask, 2)'));
      else
        assert(all(any(exon_mask([1:CFG.max_side_len-1, end-CFG.max_side_len+2:end],:),2)'));
      end  
      mask = (any(exon_mask, 2) & ~repeat_mask)';
      exon_mask = exon_mask(mask, :);
      coverage = coverage(mask, :);
    else
      error('unknown method %s', CFG.method);
    end
    if isempty(exon_mask)
      if CFG.VERBOSE>0, fprintf(1, 'no positions left for gene %i\n', g); end
      genes(g).transcript_weights(1:length(genes(g).transcripts)) = nan;
      genes(g).transcript_weights_all(1:length(genes(g).transcripts)) = nan;
      genes(g).loss.all = nan;
      genes(g).loss.exons = nan;
      genes(g).loss.introns = nan;
      continue;
    end
    
    if ~isempty(intron_mask)
      intron_count = intron_list(:,3);
    else
      intron_count = [];
    end
    
    loss_frac = zeros(length(CFG.C1_set), 1);
    weights = zeros(length(CFG.C1_set), length(gene.transcripts));
    loss = {};
    % optimise for different C1 settings
    for C1_idx = 1:length(CFG.C1_set),
      CFG.C1 = CFG.C1_set(C1_idx);
      if size(exon_mask,2)<=10
        %exon_mask = exon_mask - repmat(mean(exon_mask,1), size(exon_mask,1), 1);
        %exon_mask = exon_mask./repmat(std(exon_mask)+eps, size(exon_mask,1), 1);
        %y = full(coverage)-mean(full(coverage),1);
        
        %[weights(C1_idx,:), betas, xis, loss{end+1}] = opt_transcripts(CFG, gene, coverage, exon_mask, excluded_reads, intron_count, intron_mask, lpenv);
        %keyboard
        CFG.VERBOSE = 1;
        weights = opt_transcripts_L2(CFG, coverage, exon_mask, intron_count, intron_mask, gene.transcript_length');
      else
        weights(C1_idx,:) = nan;  
        loss{end+1} = [];
      end
    end
    [tmp, midx] = min(abs(CFG.C1_loss_frac_target-loss_frac));
    genes(g).transcript_weights = weights(midx,:);
    genes(g).transcript_weights_all = weights;
    genes(g).loss = loss{midx};
  end
end

if isequal(CFG.optimizer, 'cplex'),
  fprintf(1, '\n');
  [lpenv, status] = cplex_quit(lpenv,0);
  lpenv = 0;
end

if CFG.VERBOSE>0, fprintf(1, '\nFinished after %2.f s\n', toc); end

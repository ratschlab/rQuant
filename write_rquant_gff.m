%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% Written (W) 2009-2010 Regina Bohnert, Gunnar Raetsch
% Copyright (C) 2009-2010 Max Planck Society
%

function write_rquant_gff(genes, gff_fname, source, mapped_reads, read_len) 
% function write_rquant_gff(genes, gff_fname, source, mapped_reads, read_len) 
%
% -- input --
% genes: struct defining genes with start, stops, exons etc.
% gff_fname: name of gff file
% source: gene source
% mapped_reads: number of mapped reads
% read_len: read length  

used_contig_idx = unique([genes.chr_num]);
mapped_reads = sum(mapped_reads(used_contig_idx));
reads_per_1coverage = 1000/read_len;
rpkm_factor_all = reads_per_1coverage./(mapped_reads/10^6);
rpkm_factor_weighting = mapped_reads.*read_len;
rpkm_factor_weighting = rpkm_factor_weighting/sum(rpkm_factor_weighting);
rpkm_factor = sum(rpkm_factor_all.*rpkm_factor_weighting);

if exist(gff_fname, 'file')
  fprintf('replacing file %s...\n', gff_fname);
  [fd msg] = fopen(gff_fname, 'w+');
  disp(msg);
else
  fprintf('creating file %s...\n', gff_fname);
  [fd msg] = fopen(gff_fname, 'w+');
  disp(msg);
end

for g = 1:length(genes),
  fprintf('writing gene %i...\r', g);  
  if genes(g).is_valid==0, continue; end
  gene = genes(g);
  type  = 'gene';
  score = '.';
  phase = '.';
  for t = 1:length(gene.transcripts),
    if gene.transcript_valid(t)==0, continue; end
    type = 'transcript';
    score = '.';
    phase = '.';
    start = min(gene.exons{t}(:,1));
    stop = max(gene.exons{t}(:,2));
    % transcript
    if isfield(gene, 'transcript_weights'), % write RPKM to attribute field
      attr_str = sprintf('gene_id "%s"; transcript_id "%s"; ARC "%4.4f"; RPKM "%4.4f";', gene.name, gene.transcripts{t}, gene.transcript_weights(t), gene.transcript_weights(t)*rpkm_factor);
    else
      attr_str = sprintf('gene_id "%s"; transcript_id "%s";', gene.name, gene.transcripts{t});
    end
    fprintf(fd,'%s\t%s\t%s\t%i\t%i\t%s\t%s\t%s\t%s\n', gene.chr, source, type, start, stop, score, gene.strand, phase, attr_str);
    % exons
    attr_str = sprintf('gene_id "%s"; transcript_id "%s";', gene.name, gene.transcripts{t});
    exon_type = {'exons'};
    gff_types  = {'exon'};
    for tt = 1:length(exon_type),
      exons = gene.(exon_type{tt}){t};
      for e = 1:size(exons,1),
        type = gff_types{tt};
        score = '.';
        phase='.';
        start = exons(e,1);
        stop = exons(e,2);
        fprintf(fd,'%s\t%s\t%s\t%i\t%i\t%s\t%s\t%s\t%s\n', gene.chr, source, type, start, stop, score, gene.strand, phase, attr_str);
      end
    end
  end
end

fclose(fd);

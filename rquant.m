%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% Written (W) 2009-2010 Regina Bohnert, Gunnar Raetsch
% Copyright (C) 2009-2010 Max Planck Society
%

function rquant(CFG)
% rquant(CFG)
%
% -- input --
% CFG: configuration struct


%%%% paths
[CFG.paths CFG.rproc_par.envstr] = set_rquant_paths(CFG);

%%%% configuration
CFG = configure_rquant(CFG);

%%%% determine read length and number of mapped reads %%%%
CFG.read_len = 0;
mapped_reads = zeros(1, length(CFG.tracks_fn));
for c = 1:length(CFG.tracks_fn),
  if CFG.VERBOSE>1, fprintf(1, '\nChecking bam file: contig %i/%i\n', c, length(CFG.genome_info.flat_fnames)); end
  for f = 1:length(CFG.tracks_fn{c}),
    fname = CFG.tracks_fn{c}{f};
    assert(exist(sprintf('%s.bai', fname), 'file')==2);
    try
      [read_len tmp_mapped_reads] = get_bam_properties(CFG.tracks_fn{c}{f}, CFG.samtools_dir, CFG.genome_info.contig_names{c});
    catch
      read_len = 0;
      tmp_mapped_reads = 0;
    end
    mapped_reads(c) = mapped_reads(c) + tmp_mapped_reads;
    CFG.read_len = max(CFG.read_len, read_len);
  end
end

%%%%% rquant %%%%%
save_fname = rquant_core(CFG);

%%%%% write to gff file %%%%%
if CFG.write_gff,
  load(save_fname, 'genes');
  %unix(sprintf('cp -p %s %s/genes.mat', save_fname, CFG.out_dir)) ;
  output_file = strrep(save_fname, '.mat', 'gff3');
  write_rquant_gff(genes, output_file, 'rQuant', mapped_reads, CFG.read_len);
end

%%%%% write learned read density model %%%%%
if CFG.write_density_model,
  load(save_fname, 'RES');
  assert(length(RES)>=1);
  profiles_fn_out = strrep(save_fname, '.mat', '_profiles.txt');
  intron_dists_fn_out = strrep(save_fname, '.mat', '_intron_dists.txt');
  if CFG.max_iter > 1
    idx = length(RES)-1;
  else
    idx = 1;
  end
  write_density_model(RES{idx}.profile_weights, RES{idx}.intron_dists, profiles_fn_out, intron_dists_fn_out);
end

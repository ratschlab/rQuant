function show_genes(CFG, genes, profile_weights, intron_dists)
% show_genes(CFG, genes, profile_weights, intron_dists)
%
% -- input --
% CFG: configuration struct
% genes: struct defining a gene with start, stops, exons etc.
% profile_weights: weights of profile functions
% intron_dists: distances to closest intron

addpath('~/svn/projects/splicing/splicegraphs/');

genes = build_splice_graph_caller(genes);

for g = 1:length(genes),
  eidx = [];
  for t = 1:length(genes(g).transcripts),
    for e = 1:size(genes(g).exons{t},1),
      eidx = [eidx, genes(g).exons{t}(e,1):genes(g).exons{t}(e,2)];
    end
  end
  eidx = unique(eidx-genes(g).start+1);
  ec = zeros(genes(g).stop-genes(g).start+1,1);
  ec(eidx) = genes(g).all_coverage;
  expr = zeros(genes(g).stop-genes(g).start+1,1);
  for t = 1:length(genes(g).transcripts),
    tmp = genes(g).transcript_weights(t)*(gen_exon_features(genes(g), t, CFG.num_plifs, CFG.max_side_len)*...
                                          profile_weights(:, genes(g).transcript_len_bin(t), genes(g).expr_bin(t))) ...
                                          ;%- intron_dists(:,t);
    expr(eidx) = expr(eidx) + tmp;
  end
  [genes(g).transcript_weights sidx] = sort(genes(g).transcript_weights);
  
  genes(g).transcripts = genes(g).transcripts(sidx);
  genes(g).exons =  genes(g).exons(sidx);
  for t = 1:length(genes(g).transcripts),
    genes(g).transcripts{t} = sprintf('%10.2f    %s', genes(g).transcript_weights(t), genes(g).transcripts{t});
  end
  
  % plot splicegraph
  if isfield(genes, 'splicegraph')
    figure(); viewsplicegraph(genes(g));
    ax1 = axis;
    set(gcf, 'Position', [200 900 1500 300]);
  end
  % plot observed coverage (blue) and estimated coverage (black)
  figure(); hold on;
  plot(ec, 'b'); plot(expr, 'k'); 
  for n = 1:size(genes(g).all_introns,1),
    idx = [genes(g).all_introns(n,1):genes(g).all_introns(n,2)]-genes(g).start+1;
    plot(idx, double(genes(g).all_introns(n,3))*ones(1,length(idx)), 'r', 'LineWidth', 3);
  end
  axis tight;
  ax2 = axis;
  if isfield(genes, 'splicegraph')
    ax2(2) = ax1(2) - ax1(1); axis(ax2);
  end
  set(gcf, 'Position', [200 440 1500 400]);    
  genes(g)
  genes(g).transcript_weights
  pause
  close all;
  end
end
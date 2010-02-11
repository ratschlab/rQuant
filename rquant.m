function rquant(CFG)

%%%% paths
CFG.paths = set_rquant_paths();

%%%% configuration
CFG = configure_rquant(CFG);

if CFG.use_rproc,
  DEBUG = 0;
else
  DEBUG = 1;
  CFG.VERBOSE = 2;
end

%%%% load genes
load(CFG.gene_fn, 'genes');

% add exonic length
% initialise expression bins
% initialise transcript length bins
genes(1).expr_bin = [];
genes(1).transcript_len_bin = [];
genes(1).eidx = [];
genes(1).exonic_len = 0;
if CFG.norm_seqbias
  genes(1).num_read_starts = [];
end
del_idx = false(1, length(genes));
for g = 1:length(genes),
  eidx = [];
  min_start = inf; max_stop = -inf;
  % closed intervals
  if isequal(CFG.gene_source, 'annotation')
    genes(g).strands = repmat(genes(g).strand, 1, length(genes(g).transcripts));
  end
  if genes(g).start>genes(g).stop || genes(g).start<1 || genes(g).stop<1
    del_idx(g) = true;
  end
  genes(g).expr_bin = ones(1, length(genes(g).transcripts));
  genes(g).transcript_len_bin = ones(1, length(genes(g).transcripts));
  for t = 1:length(genes(g).transcripts),
    tidx = [];
    % closed intervals
    if isequal(CFG.gene_source, 'annotation')
      assert(isfield(genes, 'strands'));
      if genes(g).strands(t)=='+'
        genes(g).exons{t}(:,2) = genes(g).exons{t}(:,2)-1;
      else
        genes(g).exons{t}(:,1) = genes(g).exons{t}(:,1)+1;
      end
    end
    min_start = min(min(genes(g).exons{t}(:,1)), min_start);
    max_stop = max(max(genes(g).exons{t}(:,2)), max_stop);
    if any(genes(g).exons{t}(:,1)>genes(g).exons{t}(:,2)) || any(genes(g).exons{t}(:,1)<1) || any(genes(g).exons{t}(:,2)<1) || ...
      length(genes(g).exons{t}(:,1))~=length(unique(genes(g).exons{t}(:,1))) || length(genes(g).exons{t}(:,2))~=length(unique(genes(g).exons{t}(:,2))),
      del_idx(g) = true;
    end
    for e = 1:size(genes(g).exons{t},1),
      tidx = [tidx, genes(g).exons{t}(e,1):genes(g).exons{t}(e,2)];
      tidx = unique(tidx);
    end
    genes(g).transcript_length(t) = length(unique(tidx));
    if genes(g).transcript_length(t)>100000,
      del_idx(g) = true;
    end
    eidx = unique([eidx tidx]);
  end
  if isequal(CFG.gene_source, 'annotation'),
    if (min_start~=genes(g).start)
      assert(min_start==genes(g).start+1);
    end
    if (max_stop~=genes(g).stop)
      assert(max_stop==genes(g).stop-1);
    end
  end
  genes(g).start = min_start;
  genes(g).stop = max_stop;
  assert(genes(g).start<=genes(g).stop);
  genes(g).eidx = unique(eidx);
  genes(g).exonic_len = length(genes(g).eidx);
  if CFG.norm_seqbias
    genes(g).num_read_starts = zeros(1, genes(g).exonic_len);
  end
end
genes(del_idx) = [];

if CFG.subsample,
  CFG.subsample_frac_global = str2double(sprintf('%.2f', min(1,CFG.max_num_train_exm/sum([genes.exonic_len])))); 
end

% subsample here
if CFG.subsample && isequal(CFG.gene_source, 'annotation') 
  ridx = randperm(length(genes));
  CFG.subsample_idx = ridx(1:ceil(CFG.subsample_frac_global*length(genes)));
end

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
CFG.transcript_len_ranges

CFG
CFG.C2.tau
CFG.C2.kappa
CFG.C2.theta
fprintf(1,'using %i genes (deleted %i)\n\n', length(genes), sum(del_idx));
clear del_idx;

% initialise profiles and intron distances
if CFG.load_profiles,
  fprintf(1, '\nLoading profiles... (%s)\n', CFG.profiles_fn);
  load(CFG.profiles_fn, 'RES');
  profile_weights = RES{end-1}.profile_weights;
  if ~(size(profile_weights,1)==CFG.num_plifs && size(profile_weights,2)==size(CFG.transcript_len_ranges,1))
    error('profiles have wrong dimensions');
  end
  intron_dists = RES{end-1}.intron_dists;
  if ~(size(intron_dists,1)==CFG.num_intron_plifs && size(intron_dists,2)==CFG.num_intron_plifs)
    error('intron distance matrix has wrong dimensions');
  end
  clear RES;
else
  profile_weights = ones(CFG.num_plifs,1);
  x = linspace(0, 0.5, CFG.num_intron_plifs);
  intron_dists = 1-(1-x(end:-1:1))'*(1-x(end:-1:1));  
end

% iteratively optimise transcript and profile weights 
iter = 0;
while(1)
  iter = iter + 1;
  fprintf(1, '\n*** Iteration %i ***\n', iter);
  fprintf(1, '\nDetermining transcript weights...\n');
  
  PAR.CFG = CFG;
  PAR.profile_weights = profile_weights;
  PAR.intron_dists = intron_dists;
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
        %CFG.rproc_memreq = 10000;
        %CFG.rproc_par.mem_req_resubmit = [12000 20000];
       case 'drosophila',
        CFG.rproc_memreq = 4000;
        CFG.rproc_par.mem_req_resubmit = [7000 10000 20000];
        %CFG.rproc_memreq = 10000;
        %CFG.rproc_par.mem_req_resubmit = [12000 20000];
       case 'human',
        CFG.rproc_memreq = 5000;
        CFG.rproc_par.mem_req_resubmit = [8000 12000 20000];
        %CFG.rproc_memreq = 10000;
        %CFG.rproc_par.mem_req_resubmit = [12000 20000];
      end
      CFG.rproc_par.identifier = sprintf('rq.%s.f%i-', CFG.organism(1:3), f);
      fprintf(1, 'Submitting job %i (%s) to cluster\n', f, CFG.rproc_par.identifier);
      JOB_INFO(f) = rproc('opt_transcripts_caller', PAR, CFG.rproc_memreq, CFG.rproc_par, CFG.rproc_time); pause(1);
    end
    save(sprintf('~/tmp/%s_%s_%s.mat', CFG.exp, CFG.gene_source, datestr(now,'yyyy-mm-dd_HHhMM') ), 'JOB_INFO');
    % wait for jobs
    [JOB_INFO num_crashed] = rproc_wait(JOB_INFO, 60, 1, -1);
    if num_crashed>0, pause(60); end;% some delay to wait until results are written
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

  tmp = [genes.loss]; tmp = [tmp(:)];
  loss.all = [tmp.all]; loss.exons = [tmp.exons]; loss.introns = [tmp.introns];
  RES{iter}.median_loss = median(loss.all(~isnan(loss.all)));
  RES{iter}.median_loss_frac = median(loss.all(~isnan(loss.introns))./loss.exons(~isnan(loss.exons)));
  fprintf(1, '\nMedian loss after iteration %i: %5.2f\n', iter, RES{iter}.median_loss);
  fprintf(1, '\nMedian ratio loss introns/exons after iteration %i: %1.3f\n', iter, RES{iter}.median_loss_frac);
  
  if DEBUG, keyboard; end
  
  if isequal(CFG.gene_source, 'annotation')
    save_fname = sprintf('%s%s_%s_%s_iter%i.mat', CFG.out_dir, CFG.exp, CFG.gene_source, CFG.method, iter);
  else
    if CFG.paired,
      save_fname = strrep(CFG.gene_fn, '.mat', '_rquant_profile_paired.mat');
    else
      save_fname = strrep(CFG.gene_fn, '.mat', '_rquant_profile.mat');
    end
  end
  genes_tmp = genes;
  %genes = rmfield(genes, 'eidx');
  save(save_fname, 'CFG', 'RES', 'genes');
  genes = genes_tmp;
  clear genes_tmp;

  if iter >= CFG.max_iter
    fprintf(1, '\n\nFinished transcript quantitation.\n');
    break;
  end
  
  if CFG.norm_seqbias
    fprintf(1, '\nDetermining sequence bias...\n\n');
    take_idx = false(1, length(genes));
    tw = [genes.transcript_weights];
    tw = tw(~isnan(tw) | tw>0);
    min_expr = prctile(tw, 90);
    for g = 1:length(genes),
      if sum(genes(g).transcript_weights>min_expr)==1
        take_idx(g) = true;
      end
    end
    assert(sum(take_idx)>0);
    seq_norm_genes = genes(take_idx);
    for g = 1:length(seq_norm_genes),
      t = find(seq_norm_genes(g).transcript_weights>min_expr);
      assert(length(t)==1)
      seq_norm_genes(g).transcripts = seq_norm_genes(g).transcripts{t};
      seq_norm_genes(g).exons = seq_norm_genes(g).exons{t};
    end
    fprintf('using %i genes for sequence normalisation\n', sum(take_idx));
    CFG.RR.seq_norm_weights = train_norm_sequence(CFG, seq_norm_genes);
    if DEBUG, keyboard; end
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
  
  if DEBUG,
    keyboard;
    figure(); plot(profile_weights);
    figure(); imagesc(intron_dists); colorbar;
  end

  if iter>1,
    RES{iter-1}.genes = [];
    RES{iter-1}.profile_genes = [];
  end
  RES{iter}.profile_weights = profile_weights;
  RES{iter}.intron_dists = intron_dists;
  RES{iter}.profile_genes = profile_genes;
  
  save_fname = sprintf('%s%s_%s_%s_iter%i.mat', CFG.out_dir, CFG.exp, CFG.gene_source, CFG.method, iter);
  save_fname = sprintf('%s%s_%s_%s_iter%i_interm.mat', CFG.out_dir, CFG.exp, CFG.gene_source, CFG.method, iter);
  genes_tmp = genes;
  genes = rmfield(genes, 'eidx');
  RES{iter}.genes = genes;
  save(save_fname, 'CFG', 'RES', 'genes');
  genes = genes_tmp;
  clear genes_tmp;
end
addpath('~/svn/tools/rproc');
addpath('~/svn/tools/utils');
addpath('~/svn/tools/genomes');

organisms = {'elegans', 'human'};

exp{1} = {'CIseli_full', 'GRaetsch_filtered', 'GRaetsch_full', 'LPachter_full', 'MStanke_filtered', 'MStanke_full', 'SWhite_full', 'TAlioto_full', 'TWu_full'};
exp{2} = {'CIseli_full', 'GRaetsch_filtered', 'GRaetsch_full', 'LPachter_full', 'MStanke_filtered', 'MStanke_full', 'SWhite_full', 'TAlioto_full', 'TWu_full', 'MGerstein_full'};

subexp = {'sanitized.region_filtered', 'sanitized.unique.region_filtered'}; % 'sanitized.region_filtered' 'sanitized.region_filtered.clip_filtered' 'sanitized.unique.region_filtered' 'sanitized.unique.region_filtered.clip_filtered'

CFG.read_len = 75;
CFG.gene_source = 'annotation';
max_mismatches = [0:10];

%%%%% result directory %%%%
CFG.out_dir = '';

%%%%% optimizer %%%%%
CFG.optimizer = 'cplex';

%%%%% number of iterations (1: no profile learning) %%%%%
CFG.max_iter = 1;

%%%%% pre-learned profiles %%%%%
CFG.profiles_fn = '';

%%%%% rproc settings for rquant subjobs %%%%%
CFG.use_rproc = 0; % 1: cluster submission or 0: locally

%%%%% rproc settings for main job %%%%%
rproc_memreq                = 1700;
rproc_par.priority          = 8;
rproc_par.express           = 0;
rproc_par.immediately_bg    = 0;
rproc_par.immediately       = 0;
rproc_par.arch              = 64;
rproc_par.identifier        = '';
rproc_par.verbosity         = 0;
rproc_time                  = 72*60; % mins
rproc_par.envstr            = 'export MOSEKLM_LICENSE_FILE=/fml/ag-raetsch/share/software/mosek/6/licenses/mosek.lic; export LD_LIBRARY_PATH=/fml/ag-raetsch/share/software/mosek/6/tools/platform/linux64x86/bin';

cnt = 1;
for o = 1:length(organisms),
  CFG.organism = organisms{o};
  
  %%%%% directories from which to load read data and genes %%%%%
  switch CFG.organism
   case 'human'
    CFG.tracks_dir = '/fml/ag-raetsch/nobackup/projects/rgasp.2/eval/current_submissions/human_laneFiltered/';
    CFG.repeats_fn = '/fml/ag-raetsch/nobackup/projects/rgasp.2/annotations/human/hg19/repeat_masker/tracks/';
   case 'elegans'
    CFG.tracks_dir = '/fml/ag-raetsch/nobackup/projects/rgasp.2/eval/current_submissions/elegans/';
    CFG.repeats_fn = '/fml/ag-raetsch/nobackup/projects/rgasp.2/annotations/elegans/repeat_masker/tracks/';
  end
  
  %%%%% genes %%%%%
  switch CFG.organism
   case 'human'
    CFG.gene_fn = '/fml/ag-raetsch/share/projects/rquant/data_nanostring/rgasp/nanostr_genes_HepG2.mat';
   case 'elegans'
    CFG.gene_fn = '/fml/ag-raetsch/share/projects/rquant/data_nanostring/rgasp/nanostr_genes_C_elegans_L3.mat';
  end
  
  %%%%% genome info %%%%%
  switch CFG.organism
   case 'human'
    CFG.genome_info = init_genome('/fml/ag-raetsch/nobackup/projects/rgasp/genomes/human/hg19/hg19.gio/genome.config');
   case 'elegans'
    CFG.genome_info = init_genome('/fml/ag-raetsch/nobackup/projects/rgasp/genomes/elegans/elegans.gio/genome.config');
  end
  % contig length
  for c = 1:length(CFG.genome_info.flat_fnames),
    d = dir(CFG.genome_info.flat_fnames{c});
    CFG.chr_len(c) = d(1).bytes;
  end
  
  for e = 2%1:length(exp{o}),
    CFG.exp = exp{o}{e};
    for s = 1:length(subexp),
      fprintf(1, '*** %s %s %s ***\n', CFG.organism, CFG.exp, subexp{s});
      %%%%% alignments %%%%%
      for c = 1:length(CFG.genome_info.flat_fnames),
        CFG.tracks_fn{c} = {sprintf('%s%s/%s.bam', CFG.tracks_dir, CFG.exp, subexp{s})};
      end
      switch CFG.organism
       case 'human'
        CFG.tracks_max_intron_len = 200000; 
       case 'elegans'
        CFG.tracks_max_intron_len = 20000;
      end
      CFG.tracks_min_exon_len = 15;
      CFG.tracks_max_mismatches = 0;
      
      %for n = 1:length(min_exon_len),
      for n = 1:length(max_mismatches),
        CFG.tracks_max_mismatches = max_mismatches(n);
        %CFG.tracks_min_exon_len = min_exon_len(n);
        %CFG.out_fn = strrep(CFG.tracks_fn{1}{1}, '.bam', sprintf('.el_%i.mat', CFG.tracks_min_exon_len));
        CFG.out_fn = strrep(CFG.tracks_fn{1}{1}, '.bam', sprintf('.mm_%i.mat', CFG.tracks_max_mismatches));
        if 1
          rproc_par.identifier = sprintf('rq.%s%s-%i-', CFG.organism(1:2), CFG.exp(1:2), max_mismatches(n));
          fprintf(1, 'Submitting job %i (%s) to cluster\n', cnt, rproc_par.identifier);
          job(cnt) = rproc('rquant', CFG, rproc_memreq, rproc_par, rproc_time);
          cnt = cnt + 1;
        else
          rquant(CFG);
        end
      end
    end
  end
end
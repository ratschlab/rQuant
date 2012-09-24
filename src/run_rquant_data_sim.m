addpath('~/svn/tools/rproc');
addpath('~/svn/tools/utils');
addpath('~/svn/tools/genomes');

%%%%% experiment %%%%%
CFG.organism = 'human' ;% 'elegans'; % 'elegans', 'arabidopsis', 'human'
                        %CFG.exp = 'fs_weak_bias_paired'; % 'fs_strong_bias', 'fs_strong_bias_seq_bias', 'fs_weak_bias_paired'
CFG.exp = 'fs_strong_bias' %, 'fs_strong_bias_seq_bias', 'fs_weak_bias_paired'
PAR.CFG.read_len = 75;

%%%%% tracks, repeats, genes, genome info %%%%% 
switch CFG.organism
 case 'arabidopsis'
  CFG.base_dir = '/fml/ag-raetsch/share/projects/rquant/data_sim/arabidopsis/TAIR10';
  PAR.CFG.repeats_fn = '';
  PAR.CFG.genome_info = init_genome('/fml/ag-raetsch/share/databases/genomes/A_thaliana/arabidopsis_tair10/annotations/genome.gio/genome.config');
  PAR.CFG.correct_intervals = 1;
 case 'elegans'
  CFG.base_dir = '/fml/ag-raetsch/share/projects/rquant/data_sim/elegans/WS200';
  PAR.CFG.repeats_fn = '/fml/ag-raetsch/nobackup/projects/rgasp.2/annotations/elegans/repeat_masker/tracks';
  PAR.CFG.genome_info = init_genome('/fml/ag-raetsch/nobackup/projects/rgasp/genomes/elegans/elegans.gio/genome.config');
  PAR.CFG.correct_intervals = 1; % 1 if genes is half open -> corrects to closed
 case 'human'
  CFG.base_dir = '/fml/ag-raetsch/share/projects/rquant/data_sim/human/HG19';
  PAR.CFG.repeats_fn = '';
  PAR.CFG.correct_intervals = 1;
end
%PAR.anno_dir = sprintf('%s/annotation', CFG.base_dir); % annotation, name must be genes.mat
PAR.anno_dir = sprintf('%s/annotation/%s', CFG.base_dir, CFG.exp);
%PAR.track = sprintf('%s/tracks/%s.bam', CFG.base_dir, CFG.exp)
%PAR.track = sprintf('%s/tracks/%s.bam', CFG.base_dir, CFG.exp); % alignents (bai must be present)
%PAR.track = '/fml/ag-raetsch/nobackup/projects/mip_spladder/alignments/human/artifical_reads_ns/reads_strong_bias.noise0.015.sorted.bam'
%PAR.track = '/fml/ag-raetsch/nobackup/projects/mip_spladder/alignments/human/artifical_reads_ns/reads_strong_bias.sorted.bam'
%PAR.track = '/fml/ag-raetsch/nobackup/projects/mip_spladder/alignments/human/artifical_reads_ns/reads_strong_bias.noise0.015.mmr.sorted.bam' ;
%PAR.track = '/fml/ag-raetsch/nobackup/projects/mip_spladder/alignments/human/artifical_reads_ns/reads_weak_bias.noise0.015.mmr.sorted.bam' ;
%PAR.track = '/fml/ag-raetsch/nobackup/projects/mip_spladder/alignments/human/artifical_reads_ns/reads_strong_bias.noise0.015.best.sorted.bam' ;
%PAR.track = '/fml/ag-raetsch/nobackup/projects/mip_spladder/alignments/human/artifical_reads_ns/reads_weak_bias.noise0.015.sorted.bam' ;
%PAR.track = '/fml/ag-raetsch/nobackup/projects/mip_spladder/alignments/human/artifical_reads_ns/reads_strong_bias.noise0.015.cross3.bam' ;

PAR.track = '/fml/ag-raetsch/nobackup/projects/mip_spladder/alignments/human/artifical_reads_ns/reads_weak_bias.noise0.03.x2.sorted.bam' ;
%PAR.CFG.repeats_fn='/fml/ag-raetsch/nobackup/projects/mip_spladder/alignments/human/artifical_reads_ns/reads_strong_bias.noise0.015.mmr1.rep/' ;
%PAR.CFG.repeats_fn='/fml/ag-raetsch/nobackup/projects/mip_spladder/alignments/human/artifical_reads_ns/reads_strong_bias.noise0.03.x2.mmr1.rep/' ;

%%%%% output files %%%%%
PAR.output_dir = '';
PAR.output_file = '';
PAR.CFG.write_gff = 0;
PAR.CFG.write_density_model = 0;
PAR.profiles_fn_out = '';

%%%%% profile learning %%%%%
% enables profile learning
PAR.learn_profiles = 0; % 0: no learning, 1: empirically estimated, 2: optimised
% pre-learned profiles
PAR.load_profiles = 0;
PAR.profiles_fn = '';
% regularisation strengths 
C_I = 100; %[10^0 10^1 10^2 10^3];
C_F = 10^5; %[10^3 10^4 10^5 10^7]; (only used for profile optimisation)
C_N = 1; %[10^0 10^2 10^5]; (only used for profile optimisation)
C_PE = [10^2];%[10 10^2 10^3 10^4 5*10^4];

%%%%% paired-end %%%%%
PAR.CFG.paired = 0;

%%%%% sequence bias normalisation %%%%%
PAR.CFG.norm_seqbias = 0;

%%%%% rproc settings for rquant subjobs %%%%%
PAR.CFG.use_rproc = 0; % 1: cluster submission or 0: locally
if PAR.CFG.use_rproc,
  PAR.CFG.rproc_num_jobs              = 20;
  PAR.CFG.rproc_memreq                = 3000;
  PAR.CFG.rproc_par.priority          = 8;
  PAR.CFG.rproc_par.resubmit          = 3;
  PAR.CFG.rproc_par.mem_req_resubmit  = [8000 12000 20000];
  PAR.CFG.rproc_par.time_req_resubmit = [36*60 70*60 90*60];
  PAR.CFG.rproc_par.express           = 0;
  PAR.CFG.rproc_par.immediately_bg    = 0;
  PAR.CFG.rproc_par.immediately       = 0;
  PAR.CFG.rproc_par.arch              = 64;
  PAR.CFG.rproc_par.identifier        = '';
  PAR.CFG.rproc_par.verbosity         = 1;
  PAR.CFG.rproc_time                  = 5*60;
end

% !!!
PAR.CFG.tracks_max_mismatches=0

run_local = 1;

for p = 1:length(C_PE),
PAR.CFG.C_PE = C_PE(p);
for s = 1:length(C_I),
  PAR.CFG.C_I = C_I(s);
  for f = 1:length(C_F),
    PAR.CFG.C_F = C_F(f);
    for n = 1:length(C_N),
      PAR.CFG.C_N = C_N(n);
      %%%%% result directory %%%%%
      %date_exp = datestr(now,'yyyy-mm-dd');
      date_exp = datestr(now,'yyyy-mm-dd_HHhMM')
      switch CFG.organism
       case 'arabidopsis'
        PAR.output_dir = sprintf('/fml/ag-raetsch/share/projects/rquant/data_sim/arabidopsis/TAIR10/rquant/%s_%s', CFG.exp, date_exp);
       case 'elegans'
        PAR.output_dir = sprintf('/fml/ag-raetsch/share/projects/rquant/data_sim/elegans/WS200/rquant/%s_%s', CFG.exp, date_exp);
       case 'human'
        PAR.output_dir = sprintf('/fml/ag-raetsch/share/projects/rquant/data_sim/human/HG19/rquant/%s_%s', CFG.exp, date_exp);
      end
      PAR.output_file = sprintf('%s/%s_rquant.gff3', PAR.output_dir, CFG.exp);
      if ~exist(PAR.output_dir ,'dir'),
        [s m mid] = mkdir(PAR.output_dir);
        assert(s);
      end
      PAR.profiles_fn = sprintf('%s/profiles.mat', PAR.output_dir);
      % save parameters
      fname = strrep(PAR.output_file, '.gff3', '.par');
      write_parameters(PAR, fname);
      %%%%% rproc settings for main job %%%%%
      if ~run_local
        rproc_memreq                = 2000;
        rproc_par.priority          = 8;
        rproc_par.express           = 0;
        rproc_par.immediately_bg    = 0;
        rproc_par.immediately       = 0;
        rproc_par.arch              = 64;
        rproc_par.identifier        = '';
        rproc_par.verbosity         = 1;
        rproc_time                  = 72*60;
        rproc_par.identifier = sprintf('rq.%s-', CFG.organism(1:2));
        fprintf(1, 'Submitting job %s to cluster\n', rproc_par.identifier);
        job = rproc('rquant_rproc', PAR, rproc_memreq, rproc_par, rproc_time);
        pause(60);
      else
        rquant_rproc(PAR);
        out_fn=strrep(PAR.output_file, '.gff3', '.mat')
        evaluation(out_fn, 'Pearson', 2, 1) ;
      end
    end
  end
end
end


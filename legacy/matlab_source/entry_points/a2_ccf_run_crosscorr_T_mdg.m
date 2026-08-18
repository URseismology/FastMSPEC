%% This is the Entry point of the CCF code. 
%%% Declare Results Path and Environments and Component (Love/Raleigh) to Calculate the CCFs
clear;clc;

env = 'test';
comp = 'love/madagascar';
results_root_path = ['/scratch/tolugboj_lab/Sayan_Swar_WS/PythonEnv/Python_Notebooks/PRJ_SPAC/results/' env '/' comp '/'];
data_root_path = ['/scratch/tolugboj_lab/Sayan_Swar_WS/PythonEnv/Python_Notebooks/PRJ_SPAC/data/' env '/'];

%% Declare Experimental Setup
IsFigure1 = 1;
IsFigure2 = 0;

% OUTPUT SETTINGS
stacktype.IsOutputFullstack = 1; %save full year ccf stacks
stacktype.IsOutputMonthstack = 0; %save month ccf stacks
stacktype.IsOutputDaystack = 0; %save day ccf stacks
stacktype.IsOutputSinglestack = 0; %save single ccf before stacking
stacktype.IsOutputSeismograms = 0; %save raw seismograms before cross-correlating

% GENERAL PROCESSING
filttype.dt=1;
filttype.winlength = 3; %segmenting each day record in N hours with overlap
filttype.IsRemoveIR = 0; %remove instrument response
filttype.units_RemoveIR = 'M'; %'M' displacement | 'M/S' velocity
filttype.IsDetrend = 0; %detrend the data
filttype.IsTaper = 0; %Apply cosine taper to data chunks

% PERFORMANCE EVALUATION
filttype.Iseval = 0; %Make 1 if running for performance evaluation and increase window length
filttype.Chunklen = 300;
zero_padding_list = [8192,28800,57600,86400,115200,129600,172800];

% ADVANCED PREPROCESSING
% (1) ONE-BIT NORMALIZATION & SPECTRAL WHITENING? (Bensen et al. 2007)
filttype.IsSpecWhiten = 0; %Whiten spectrum
filttype.IsOBN = 0; %One-bit normalization

% (2) TIME-FREQUENCY NORMALIZATION (Ekstrom et al. 2009; Shen et al. 2011)
filttype.IsFTN = 0; %frequency-time normalization? (If 1, applied instead of whitening and one-bit normalization)
filttype.frange_FTN = [1/60 1/10]; %frequency range over which to construct FTN seismograms

filttype.IsMultiTaper = 0; %old code with Mspec, same is available with filttype variable

% (3) BASIC PREFILTER (Ekstrom 2011)
filttype.IsPrefilter = 0; %apply butterworth bandpass filter before cross-correlation?
filttype.frange_prefilt = [1/100 1/10];

% (4) CODA PROCESSING (https://doi.org/10.1016/j.earscirev.2020.103285)
filttype.CodaProcessing = 0;

% (4) MULTITAPER SPECTRAL CALCULATIONS 
filttype.IsMspec = 1; %new version
filttype.technique = 'FastMspec'; %choose technique: FastMspec, Mspec, MspecBestK

filttype.Wband = 0.001; %0.002; %2.0e-3;0.01;
filttype.epsilon = 1e-5;
filttype.cutoff = 1-filttype.epsilon;
%filttype.cutoff = 1-1e-5; %1-1e-9; %1-1e-16; %

if (filttype.IsMspec == 1) & (filttype.technique == 'FastMspec')
    disp('Computing Multi Tapers')
    FMTSE = FastMultitaper(filttype.dt*filttype.winlength*60*60+1,filttype.Wband,filttype.cutoff,filttype.epsilon);
    filttype.FMTSE = FMTSE;
end

filttype.NW_mspec = 100; %Time-Bandwidth Product
filttype.K_taps_mspec = 80; %2*filttype.NW - 1; %Number of Tapers

expname = 'fastmspec';

%% Load File List for CCF Computation and Setup Results Directory
ccflist_filepath = '/scratch/tolugboj_lab/Sayan_Swar_WS/PythonEnv/Python_Notebooks/PRJ_SPAC/data/test/metadata/madagascar_stn_conn_ccflist.csv';
ccflist_stnconn = readtable(ccflist_filepath);
ccflist_stnconn = ccflist_stnconn(ccflist_stnconn.filesize_mb>100,:);
ccflist_filenames = ccflist_stnconn.filelocation;
nconn = size(ccflist_stnconn,1);

parameters = ccf_setup_params_T_mdg(data_root_path,comp,filttype); 
path_vars = ccf_setup_results_directory_mdg(results_root_path,parameters,stacktype,expname);
%ccf_ready_data = load(ccflist_filenames{1});

%%
code_start2 = tic;
IsFigure1 = 1;

for i=1:nconn
    filenm = ccflist_filenames{i};
    ccf_ready_data = load(filenm);

    ccf_ready_data = ccf_preprocess_filter_data(filttype,ccf_ready_data,parameters); %perform filtering/processing on the frequency domain
    saved_ccf_path = ccf_compute_crosscorr_T(ccf_ready_data,path_vars,parameters,stacktype,filttype); %calculate cross correlation
    
    if IsFigure1
        smooth_val = 1; %5
        ccf_basic_plot_time_freq(saved_ccf_path.fullstack,parameters,ccf_ready_data.stapairsinfo,path_vars,smooth_val,1); %plot results single station ccf experiment 
    end

    ccf_compute_details(i).saved_ccf_path = saved_ccf_path.fullstack;
    ccf_compute_details(i).stapairsinfo = ccf_ready_data.stapairsinfo;
    
    %if i==2
    %   break
    %end

end
code_end2 = toc(code_start2);

ccfcompute_file_name = [path_vars.exp_dir,'/ccf_compute_log.mat'];
save(ccfcompute_file_name, 'ccf_compute_details');

total_time = code_end2;
ccf_save_paramdetails(path_vars,parameters,filttype,stacktype,total_time); %save all parameters for the experiment performed


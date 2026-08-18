function [parameters] = ccf_setup_params_Z(data_root_path,comp,filttype)

    addpath('/scratch/tolugboj_lab/Sayan_Swar_WS/PythonEnv/Python_Notebooks/PRJ_SPAC/codes/test/matlab/functions/');
    addpath('/scratch/tolugboj_lab/Sayan_Swar_WS/PythonEnv/Python_Notebooks/PRJ_SPAC/codes/test/matlab/functions/calc_Rayleigh_disp/');
    addpath('/scratch/tolugboj_lab/Sayan_Swar_WS/PythonEnv/Python_Notebooks/PRJ_SPAC/codes/test/matlab/functions/jCommon/');
    addpath('/scratch/tolugboj_lab/Sayan_Swar_WS/PythonEnv/Python_Notebooks/PRJ_SPAC/codes/test/matlab/functions/jSpectral/');
    addpath('/scratch/tolugboj_lab/Sayan_Swar_WS/PythonEnv/Python_Notebooks/PRJ_SPAC/codes/test/matlab/functions/jVarfun/');
    
    %%% --- Data/Metadata Path Declarations --- %%%
    [stalist, stalat, stalon, staz] = textread([data_root_path 'metadata/sta_list.txt'],'%s %f %f %f\n');
    parameters.stalist = stalist;
    parameters.stalat = stalat;
    parameters.stalon = stalon;
    parameters.staz = staz;
    parameters.nsta = length(parameters.stalist);    
    parameters.PZpath = '../INSTRUMENT/'; % path to RESP files containing poles and zeros
    parameters.orientation_path = './OBS_orientations.txt'; % Column 1: station name;   Column 2: H1 degrees CW from N
    
    parameters.datapath = [data_root_path 'raw_data/']; 
    parameters.proc_datapath = [data_root_path 'processed_data/' comp '/']; 

    A = readtable([data_root_path 'metadata/orientation.csv']);
    parameters.slist = A{:,'Var1'};
    parameters.orientations = A{:,'Var2'};

    clear stalist stalat stalon staz
    
    %%% --- Parameters to build up gaussian filters --- %%%
    parameters.min_width = 0.18;
    parameters.max_width = 0.30;
    
    %%% --- Parameters for initial processing --- %%%
    parameters.dt = filttype.dt; % sample rate
    parameters.comp = 'NA'; % component %not being used in any code
    parameters.mindist = 20; % min. distance in kilometers
    parameters.year = '';

    %%% --- Parameters for ccf_ambnoise --- %%%
    parameters.winlength = filttype.winlength; %hours
    parameters.Nstart_sec = 50; % number of sections to offset start of seismogram
    parameters.Nstart = parameters.Nstart_sec/parameters.dt;
    
    %%% --- Parameters for fitbessel --- %%%
    parameters.npts = parameters.winlength*3600 / parameters.dt;
    
    %%% --- Parameters for using Radon Transform picks --- %%%
    parameters.path_LRT_picks = './mat-LRTdisp/LRT_picks/';
    
    %%% --- Parameters to Define Channel and Braodband Data Type --- %%%
    parameters.strSACcomp = 'Z';
    parameters.strNAMEcomp = 'ZZ';

    parameters.strSACcomp1 = 'N';
    parameters.strNAMEcompR = 'RR';
    
    parameters.strSACcomp2 = 'E';
    parameters.strNAMEcompT = 'TT';
    
    %%% --- Parfor Setup --- %%%
    parameters.Nworkers = 5;
    parameters.comp=comp;




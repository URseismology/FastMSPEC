%%  Script to parallelize computing CCF
%   Authors: Siyu Xue


% connect to the cluster
d =  parcluster('Bluehive_r2019a'); 
d.AdditionalProperties.AdditionalSubmitArgs='-t 7200 -p urseismo'; %-t in minutes
d.NumWorkers = 250; %set number of workers
ap = {'/scratch/tolugboj_lab/Prj2_SEUS_RF/3_Src'};


% ---- read in the station pair list ------
StationList = '/gpfs/fs2/scratch/tolugboj_lab/Prj5_HarnomicRFTraces/Extra_from_noise/CCF_auto/connected_pairs.txt';
%[stalist1, stanet1, stacomp1, stalat1, stalon1, staz1, stalist2, stanet2, stacomp2, stalat2, stalon2, staz2] = textread([StationList],'%s %s %s %f %f %f %s %s %s %f %f %f\n'); 
[stanet1, stalist1, stacomp1, stanet2, stalist2, stacomp2] = textread([StationList],'%s %s %s %s %s %s\n'); 
nsta = length(stalist1);
winlength = 4;

% ------ set some paths ------
parameters.workingdir = '/gpfs/fs2/scratch/tolugboj_lab/Prj5_HarnomicRFTraces/Extra_from_noise/CCF_auto/';

parameters.ccfpath = [parameters.workingdir,'ccf/'];
parameters.figpath = [parameters.workingdir,'figs/'];
parameters.seis_path = [parameters.workingdir,'seismograms/'];


% ------ Build File Structure: cross-correlations -------
ccf_path = parameters.ccfpath;
ccf_winlength_path = [ccf_path,'window',num2str(winlength),'hr/'];
ccf_singlestack_path = [ccf_winlength_path,'single/'];
ccf_daystack_path = [ccf_winlength_path,'dayStack/'];
ccf_monthstack_path = [ccf_winlength_path,'monthStack/'];
ccf_fullstack_path = [ccf_winlength_path,'fullStack/'];

if ~exist(ccf_path)
    mkdir(ccf_path)
end
if ~exist(ccf_winlength_path)
    mkdir(ccf_winlength_path)
end
if ~exist(ccf_singlestack_path)
    mkdir(ccf_singlestack_path)
end
if ~exist(ccf_daystack_path)
    mkdir(ccf_daystack_path)
end
if ~exist(ccf_monthstack_path)
    mkdir(ccf_monthstack_path)
end
if ~exist(ccf_fullstack_path)
    mkdir(ccf_fullstack_path)
end

PATHS = {ccf_singlestack_path; ccf_daystack_path; ccf_monthstack_path; ccf_fullstack_path};
for ipath = 1:length(PATHS)
    ccfR_path = [PATHS{ipath},'ccfRR/'];
    ccfT_path = [PATHS{ipath},'ccfTT/'];
    ccfZ_path = [PATHS{ipath},'ccfZZ/'];
    if ~exist(ccfR_path)
        mkdir(ccfR_path);
    end
    if ~exist(ccfT_path)
        mkdir(ccfT_path);
    end
    if ~exist(ccfZ_path)
        mkdir(ccfZ_path);
    end
end

% Build File Structure: figures
figpath = parameters.figpath;
fig_winlength_path = [figpath,'window',num2str(winlength),'hr/'];
if ~exist(figpath)
    mkdir(figpath);
end
if ~exist(fig_winlength_path)
    mkdir(fig_winlength_path);
end

% Build File Structure: windowed seismograms
seis_path = parameters.seis_path;
seis_winlength_path = [seis_path,'window',num2str(winlength),'hr/'];
if ~exist(seis_path)
    mkdir(seis_path);
end
if ~exist(seis_winlength_path)
    mkdir(seis_winlength_path);
end
parameters.figpath = [parameters.workingdir,'figs/'];


% ------ parallel the job ------

for istai = 1:nsta % parallel each station pair
    allJobs = [];
    sta1=char(stalist1(istai,:));
    comp1=char(stacomp1(istai,:)); 
    net1=char(stanet1(istai,:)); 
    sta2=char(stalist2(istai,:));
    comp2=char(stacomp2(istai,:)); 
    net2=char(stanet2(istai,:)); 

    remoteJob = batch(d, @a1_ccf_ambnoise_RTZ_NE_Para, 0, {sta1, comp1, net1, sta2, comp2, net2, PATHS}, ...
    'AdditionalPaths', ap, ...
    'CaptureDiary', false, 'AutoAddClientPath', false, 'Pool', 2); % get each job 2 workers

    allJobs = [allJobs remoteJob];
        
 
end
 

  function [path_vars] = ccf_setup_results_directory(root_path,parameters,stacktype)
    d = dir([root_path,'exp*']);
    if length(d)>0
        max_dir = max(cellfun(@(x) str2double(regexp(x, '\d+', 'match', 'once')), {d.name}));
    else
        max_dir = 0;
    end
    
    if max_dir == 0
        dir_index = 0;
    else
        dir_index = max_dir;
    end
    
    exp_dir = [root_path, 'exp_' num2str(dir_index+1)];
    if ~exist(exp_dir, 'dir')
        mkdir(exp_dir);
    end

    ccf_path = [exp_dir, '/ccf/'];
    figpath = [exp_dir,'/figs/'];
    seis_path = [exp_dir,'/seismograms/'];   

    ccf_winlength_path = [ccf_path,'window',num2str(parameters.winlength),'hr/'];
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
    
    if stacktype.IsOutputSinglestack
        if ~exist(ccf_singlestack_path)
            mkdir(ccf_singlestack_path)
        end
        ccfZ_path = [ccf_singlestack_path,'ccf',parameters.strNAMEcomp,'/'];
        if ~exist(ccfZ_path)
            mkdir(ccfZ_path);
        end
        
    end
    
    if stacktype.IsOutputDaystack
        if ~exist(ccf_daystack_path)
            mkdir(ccf_daystack_path)
        end
        ccfZ_path = [ccf_daystack_path,'ccf',parameters.strNAMEcomp,'/'];
        if ~exist(ccfZ_path)
            mkdir(ccfZ_path);
        end
    end
    
    if stacktype.IsOutputMonthstack
        if ~exist(ccf_monthstack_path)
            mkdir(ccf_monthstack_path)
        end
        ccfZ_path = [ccf_monthstack_path,'ccf',parameters.strNAMEcomp,'/'];
        if ~exist(ccfZ_path)
            mkdir(ccfZ_path);
        end
    end
    
    if stacktype.IsOutputFullstack
        if ~exist(ccf_fullstack_path)
            mkdir(ccf_fullstack_path)
        end
        ccfZ_path = [ccf_fullstack_path,'ccf',parameters.strNAMEcomp,'/'];
        if ~exist(ccfZ_path)
            mkdir(ccfZ_path);
        end
    end
    
    
    PATHS = {ccf_singlestack_path; ccf_daystack_path; ccf_monthstack_path; ccf_fullstack_path};
    
    % Build File Structure: figures
    fig_winlength_path = [figpath,'window',num2str(parameters.winlength),'hr/'];
    if ~exist(figpath)
        mkdir(figpath);
    end
    if ~exist(fig_winlength_path)
        mkdir(fig_winlength_path)
    end
    
    % Build File Structure: windowed seismograms
    seis_winlength_path = [seis_path,'window',num2str(parameters.winlength),'hr/'];
    if ~exist(seis_path)
        mkdir(seis_path);
    end
    if ~exist(seis_winlength_path)
        mkdir(seis_winlength_path)
    end

    % Prepare Folders to Save the Final Computed CCFs
    for ista1=1:parameters.nsta
        sta1=char(parameters.stalist(ista1,:));
        for ipath = 1:length(PATHS)
            ccfZ_path = [PATHS{ipath},'ccf',parameters.strNAMEcomp,'/'];
            if ~exist([ccfZ_path,sta1])
                mkdir([ccfZ_path,sta1]);
            end

            ccfT_path = [PATHS{ipath},'ccf',parameters.strNAMEcompT,'/'];
            if ~exist([ccfT_path,sta1])
                mkdir([ccfT_path,sta1]);
            end
        end
    
        seisZ_path = [seis_winlength_path,parameters.strNAMEcomp(1),'/'];
        if ~exist([seisZ_path,sta1])
            mkdir([seisZ_path,sta1]);
        end

        seisT_path = [seis_winlength_path,parameters.strNAMEcompT(1),'/'];
        if ~exist([seisT_path,sta1])
            mkdir([seisT_path,sta1]);
        end

    end
    
    path_vars.ccf_singlestack_path = ccf_singlestack_path;
    path_vars.ccf_daystack_path = ccf_daystack_path;
    path_vars.ccf_monthstack_path = ccf_monthstack_path;
    path_vars.ccf_fullstack_path = ccf_fullstack_path;
    path_vars.fig_winlength_path = fig_winlength_path;
    path_vars.seis_winlength_path = seis_winlength_path;
    path_vars.exp_dir = exp_dir;
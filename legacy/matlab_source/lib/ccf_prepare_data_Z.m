function [saved_data,code_end] = ccf_prepare_data_Z(parameters,filttype)
    
    code_start = tic;
    stalist = parameters.stalist;
    nsta=parameters.nsta;
    dt = parameters.dt;
    winlength = parameters.winlength;
    Nstart_sec = parameters.Nstart_sec; 
    Nstart = Nstart_sec/dt;
    dist_min = parameters.mindist;
    datadir = parameters.datapath;
    strSACcomp = parameters.strSACcomp;
    proc_datapath = parameters.proc_datapath;
    year = '';
    saved_data = struct();
    PZpath = parameters.PZpath;
    total_iter=1;

    for ista1=1:nsta
    
        sta1=char(stalist(ista1,:));   
        list1 = dir([datadir,sta1,'/*',strSACcomp,'.sac']);
    
        for ista2=ista1:nsta
           
            clear lat1 lat2 lon1 lon2 dist az
            sta2=char(stalist(ista2,:));
    
            %% All Santity Testing on the Data
            % If same station, skip
            if(strcmp(sta1,sta2))
                continue
            end
            
            % If file already prepared, skip
            if exist([proc_datapath sta1 '_' sta2 '_win_' num2str(winlength) '_all_matched_data.mat'])  || ...
                exist([proc_datapath sta2 '_' sta1 '_win_' num2str(winlength) '_all_matched_data.mat'])
                disp('Data is already Prepared for this Station Pair');
                saved_data(total_iter).filename = [proc_datapath sta1 '_' sta2 '_win_' num2str(winlength) '_all_matched_data.mat'];
                total_iter = total_iter+1;
                code_end = toc(code_start);
                continue
            end
    
            % Get a list of all available data
            row_count=1;

            for ifil = 1:length(list1)
                file1cZ = list1(ifil).name;
    
                % Check that day file exists for station 2
                Nchar = length(sta1);
                file2cZ = dir([datadir,sta2,'/',sta2,file1cZ(Nchar+1:end)]);
                str = strsplit(file1cZ,'.');
                hdayid = [str{2},'.',str{3},'.',str{4},'.',str{5},'.',str{6}];
                
                if isempty(file2cZ)
                    disp(['No data for ',sta2,' on day ',hdayid,'... skipping'])
                    continue
                end
   
                clear data1cZ data2cZ
    
                disp(['Looking at ',hdayid,' ',sta2]);
                
                data1cZ= dir([datadir,sta1,'/',year,'/',sta1,'.',hdayid,'.*',strSACcomp,'.sac']);
                data2cZ= dir([datadir,sta2,'/',year,'/',sta2,'.',hdayid,'.*',strSACcomp,'.sac']);
                
                data1cZ =  [datadir,sta1,'/',year,'/',data1cZ.name]; 
                data2cZ =  [datadir,sta2,'/',year,'/',data2cZ.name]; 

  
                %------------------- TEST IF DATA EXIST------------------------
                [S1Zt,S1Zraw,S1,S1Ztstart] = load_sac(data1cZ);
                [S2Zt,S2Zraw,S2,S2Ztstart] = load_sac(data2cZ);
                
                % Check that sample rates are the same
                if S1.DELTA ~= S2.DELTA
                    error('S1 and S2 sample rates don''t match!');
                end
                
                % Make sure all times (of both waveforms) are relative to same
                % reference point. That is both waveforms start at the same time.
                starttime = S1Ztstart;
                S1Zt = S1Zt + seconds(S1Ztstart-starttime);
                S2Zt = S2Zt + seconds(S2Ztstart-starttime);
                
                % Ensure that files have same start time to within 1 sample
                if abs(seconds(S1Ztstart-S2Ztstart)) > S1.DELTA
                    error('Station files do not have same start time');
                end
                
                % Make sure sample rates all match
                if (abs(S1.DELTA-dt) >= 0.01*dt ) || (abs(S2.DELTA-dt) >= 0.01*dt )
                    error('sampling interval does not match data! check dt');
                end

                if filttype.IsRemoveIR
                    [S1Zraw, S2Zraw] = ccf_remove_instrument_response_Z(S1Zt,S2Zt,S1Zraw,S2Zraw,PZpath,sta1,sta2,units_RemoveIR='M');
                end
       
                % Check to make sure there's actual data
                zeroind1 = find(S1Zraw == 0);
                zeroind2 = find(S2Zraw == 0);
                if length(zeroind1) == length(S1Zraw) || length(zeroind2) == length(S2Zraw)
                    disp('All zeros!');
                    continue
                end
        
                if(length(S1Zt)*length(S2Zt)==0)
                    display(['no data for ! station ',sta2]);
                    continue
                end

        
                % Determine the time span to cut to ... this will change with
                % different segments
                clear tcut
        
                if length(S1Zraw) < 30000 
                    disp(['Sta1 ',sta1,' : ',num2str(length(S1Zraw)),' is too short!'])
                    continue
                elseif length(S2Zraw) < 30000 
                    disp(['Sta2 ',sta2,' : ',num2str(length(S2Zraw)),' is too short!'])
                    continue
                end
            %% All Santity Testing Ends Here. Data Preparation Begins
    
                if(~exist('lat2','var'));
        
                    lat1=S1.STLA;
                    lon1=S1.STLO;
                    dep1=S1.STEL; % depth is negative for OBS and positive for land stations
        
                    lat2=S2.STLA;
                    lon2=S2.STLO;
                    dep2=S2.STEL; % depth is negative for OBS and positive for land stations
        
                    % Get the interstation distance and azimuth
                    [delta,S1az] = distance(lat1,lon1,lat2,lon2);
                    [delta,S2az] = distance(lat2,lon2,lat1,lon1);
        
                    dist=deg2km(delta);
        
                    if(dist < dist_min)
                        display(['distance shorter than ',num2str(dist_min),' km, skip']);
                        break
                    end
                end 
    
                stapairsinfo.stanames = {sta1,sta2};
                stapairsinfo.lats = [lat1,lat2];
                stapairsinfo.lons = [lon1,lon2];
                stapairsinfo.dt = dt;
                stapairsinfo.r = dist;
               
           
                hour_length = winlength;        
                nwin = floor(24/hour_length)*2-1; 
                win_length = hour_length*60*60/dt;
                last_pt = win_length*.5*(nwin-1)+1+Nstart+win_length;
                if last_pt < length(S1Zraw)
                    nwin = nwin + 1;
                end
			    
                for iwin = 1:nwin 
                    disp([num2str(iwin),'-',hdayid]);
                    if hour_length == 24
                        pts_begin = Nstart;
                        pts_end = length(S1Zraw)-Nstart;
                    else
                        pts_begin = win_length*.5*(iwin-1)+1+Nstart;
                        pts_end = pts_begin+win_length;
                    end
        
                    if pts_begin > length(S1Zraw) || pts_begin > length(S2Zraw) || pts_end > length(S1Zraw) || pts_end > length(S2Zraw)
				        pts_begin = length(S1Zraw)-win_length-Nstart;
                        pts_end = pts_begin+win_length;
                    end
                    tcut = [pts_begin:pts_end] * dt;
        
                    % cut in time Z for STA1
                    S1Z=interp1(S1Zt,S1Zraw,tcut);
                    S1Z(isnan(S1Z))=0;
        
                    % cut in time Z for STA2
                    S2Z=interp1(S2Zt,S2Zraw,tcut);
                    S2Z(isnan(S2Z))=0;

                    S1_name = sta1; S2_name = sta2;
                    S1_data_mat(row_count,iwin,:) = S1Z;
                    S2_data_mat(row_count,iwin,:) = S2Z;
    
                end
                row_count = row_count + 1;
            end
            code_end = toc(code_start);
            matched_data_file_name = [proc_datapath sta1 '_' sta2 '_win_' num2str(winlength) '_all_matched_data.mat'];
            save(matched_data_file_name,'S1_name','S1_data_mat','S2_name','S2_data_mat','stapairsinfo')
            
            saved_data(total_iter).filename = matched_data_file_name;

            clear S1_name S1_data_mat S2_name S2_data_mat stapairsinfo matched_data_file_name
        end
        total_iter = total_iter+1;
    end


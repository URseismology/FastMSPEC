function [saved_ccf_path] = ccf_compute_crosscorr_mtc_T(ccf_data_file,path_vars,parameters,stacktype,filttype)
   
    S1_size = size(ccf_data_file.S1_data_mat);
    S2_size = size(ccf_data_file.S2_data_mat);
    stapairsinfo = ccf_data_file.stapairsinfo;
    N1 = S1_size(3); N2 = S2_size(3);

    S1_data_mat_2D = reshape(ccf_data_file.S1_data_mat, [], N1);
    S2_data_mat_2D = reshape(ccf_data_file.S2_data_mat, [], N2);

    if filttype.Iseval==1
        S1_data_mat_2D = S1_data_mat_2D(1:filttype.Chunklen,:);
        S2_data_mat_2D = S2_data_mat_2D(1:filttype.Chunklen,:);
        S1_data_mat_2D = [S1_data_mat_2D, zeros(size(S1_data_mat_2D,1),filttype.zero_padding-size(S1_data_mat_2D,2))];
        S2_data_mat_2D = [S2_data_mat_2D, zeros(size(S2_data_mat_2D,1),filttype.zero_padding-size(S2_data_mat_2D,2))];
    end
    
    clear N1 N2
    N1 = size(S1_data_mat_2D,2);
    dt=filttype.dt;

    switch filttype.technique
        case 'Mspec'
            %%%Traditional Mspec Calulcations
            %[PSI,~] = sleptap(N1, filttype.Wband*N1, ceil(2*filttype.Wband*N1) - 1);
            [PSI,~] = sleptap(N1, filttype.NW_mspec, filttype.K_taps_mspec);
            code_start=tic;
            [~,SXX,SYY,SXY,totalMB] = mspec_fast([],dt,S1_data_mat_2D',S2_data_mat_2D',PSI);
            code_end=toc(code_start);
            saved_ccf_path.technique_name = [filttype.technique '_log'];
            saved_ccf_path.taper_size = size(PSI,2);
        
        case 'FastMspec'            
            %%%(A) Thomposon Revisited Mspec Calulcations: Karnik Fused PSI           
            %FMTSE1 = FastMultitaper(N1,filttype.Wband,filttype.cutoff,filttype.epsilon);
            FMTSE1 = filttype.FMTSE;
            code_start=tic;
            [~,SXX,SYY,SXY,totalMB] = mspec_fast(FMTSE1,dt,S2_data_mat_2D',S1_data_mat_2D',FMTSE1.S);
            code_end=toc(code_start);
            saved_ccf_path.technique_name = [filttype.technique '_log'];
            saved_ccf_path.taper_size = size(FMTSE1.S,2);
        
        case 'MspecBestK'
            %%%(B) Thomposon Revisited Mspec Calulcations: Mspec PSI with Karnik Tapers
            FMTSE1 = FastMultitaper(N1,filttype.Wband,filttype.cutoff,filttype.epsilon);
            [PSI,~] = sleptap(N1,ceil(N1*filttype.Wband),FMTSE1.K);
            code_start=tic;
            [~,SXX,SYY,SXY] = mspec_fast([],dt,S1_data_mat_2D',S2_data_mat_2D',PSI);
            code_end=toc(code_start);
            saved_ccf_path.technique_name = [filttype.technique '_log'];
            saved_ccf_path.taper_size = size(PSI,2);
    end
    
    % matrix ops
    coh_sum_Z =  SXY ./ ( sqrt(SXX) .* sqrt(SYY) ); 
    coh_sum_Z(isnan(coh_sum_Z)) = 0;
    coh_sum = sum(coh_sum_Z, 2);
    coh_num_size = size(coh_sum_Z);
    coh_num = coh_num_size(2);
 
    % mspec only generative positive spectra. replicate for negative spectra
    sig_len = length(coh_sum);
    if mod(sig_len, 2) == 0
        coh_sum_neg = conj(flipud(coh_sum(2:end-1)));
    else
        coh_sum_neg = conj(flipud(coh_sum(2:end)));
    end
    coh_sum = [coh_sum; coh_sum_neg];

    sta1 = stapairsinfo.stanames{1};
    sta2 = stapairsinfo.stanames{2};
    ccfZ_fullstack_path = [path_vars.ccf_fullstack_path,'ccf',parameters.strNAMEcompT,'/'];
    save(sprintf('%s%s/%s_%s_f.mat',ccfZ_fullstack_path,sta1,sta1,sta2),'coh_sum','coh_num','stapairsinfo');
    saved_ccf_path.fullstack = sprintf('%s%s/%s_%s_f.mat',ccfZ_fullstack_path,sta1,sta1,sta2);
    if exist('K_taps_mspec','var') & exist('NW_mspec','var')
        saved_ccf_path.K_taps_mspec = K_taps_mspec;
        saved_ccf_path.NW_mspec = NW_mspec;
    end
    saved_ccf_path.runtime = code_end;
    saved_ccf_path.psi_memory_space = totalMB;

    disp(['Full Stacked Cross Correlation has been computed for ' sta1 '-' sta2  ' sation pairs'])

    clear nPos coh_sum_neg

end 


%% Additional Code for furture ref if needed
   % build full complex trace
   % Assuming coh_trace is a 1D array with length 5401 (one-sided FFT)
   %coh_trace = sum(coh_complex, 2);
   %coh_trace_full = zeros(1, S1_size(3));  % Initialize full FFT (length 10801)

   % Copy the positive frequencies from the one-sided FFT
   %coh_trace_full(1:nF) = coh_trace;  % First half is the same

   % Add the negative frequencies (conjugate symmetry)
   %coh_trace_full(nF+1:end) = conj(flip(coh_trace(2:end)));  % Flip and conjugate the positive part


   %norm_coh = coh_trace ./ max(abs(coh_trace));
   %plot(FF, norm_coh, 'linewidth', 2); title(num2str(ii))
    
    %saved_ccf_path = ccf_save_computed_ccf_Z(norm_coh,coh_num_size,coh_num,path_vars,parameters,stapairsinfo,stacktype);

    %clear coh_trace S2_data_mat S1_data_mat fftS1Z fftS2Z coh_num_size coh_num

    %plotting in case extra
        %%stack and plot just normalized real
    %fnyq= 1/2*dt;
    %df = 1/(S1_size(3));
    %FF = 0:df:fnyq; nF = length(FF);

    %sumcoh = sum(real(coh_complex), 2);
    %maxcoh = max(abs(sumcoh));
    %norm_realcoh = sumcoh./ maxcoh;
    %plot(FF, norm_realcoh, 'linewidth', 2);
function [saved_ccf_path] = ccf_compute_crosscorr_Z(ccf_data_file,path_vars,parameters,stacktype,filttype)
    
    S1_data_mat = ccf_data_file.S1_data_mat;
    S2_data_mat = ccf_data_file.S2_data_mat;
    stapairsinfo = ccf_data_file.stapairsinfo;
    
    if filttype.IsMultiTaper || filttype.IsFTN || filttype.IsOBN || filttype.IsSpecWhiten
        fftS1Z = S1_data_mat;
        fftS2Z = S2_data_mat;
    
    elseif filttype.IsMspec
        saved_ccf_path = ccf_compute_crosscorr_mtc_Z(ccf_data_file,path_vars,parameters,stacktype,filttype);
        return
    
    else
        fftS1Z = fft(S1_data_mat,[],3);
        fftS2Z = fft(S2_data_mat,[],3);
    end

    coh_trace = fftS2Z .* conj(fftS1Z);
    coh_num_size = size(coh_trace);
    coh_num = coh_num_size(1)*coh_num_size(2);
    coh_trace = coh_trace ./ abs(fftS1Z) ./ abs(fftS2Z);
    coh_trace(isnan(coh_trace)) = 0;
    
    saved_ccf_path = ccf_save_computed_ccf_Z(coh_trace,coh_num_size,coh_num,path_vars,parameters,stapairsinfo,stacktype);

    clear coh_trace S2_data_mat S1_data_mat fftS1Z fftS2Z coh_num_size coh_num

end 


%% Extra Code for furture ref if needed
%%% Computing CCF for Sta2 - Sta1, Swaping the stations
%sta1 = stapairsinfo{2};
%sta2 = stapairsinfo{1};

%fftS1Z = fft(S2_data_mat,[],3);
%fftS2Z = fft(S1_data_mat,[],3);

%coh_trace = fftS2Z .* conj(fftS1Z);
%coh_num = length(coh_trace);
%coh_trace = coh_trace ./ abs(fftS1Z) ./ abs(fftS2Z);
%coh_trace(isnan(coh_trace)) = 0;

%coh_trace_sum = sum(sum(coh_trace,1),2);
%coh_sumZ = reshape(coh_trace_sum,1,[]);

%ccfZ_fullstack_path = [path_vars.ccf_fullstack_path,'ccf',parameters.strNAMEcomp,'/'];
%save(sprintf('%s%s/%s_%s_f.mat',ccfZ_fullstack_path,sta1,sta1,sta2),'coh_sumZ','coh_num','stapairsinfo');
    
    
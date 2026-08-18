function [ccf_ready_data] = ccf_preprocess_filter_data(filttype,ccf_ready_data,parameters)
    
    S1_size = size(ccf_ready_data.S1_data_mat);
    S2_size = size(ccf_ready_data.S2_data_mat);

   if filttype.IsDetrend
        ccf_ready_data.S1_data_mat = ccf_detrend_3dim(ccf_ready_data.S1_data_mat,S1_size);        
        ccf_ready_data.S2_data_mat = ccf_detrend_3dim(ccf_ready_data.S2_data_mat,S2_size);    
        disp('Detrending operation on the signal completed')
   end

   if filttype.IsTaper
        ccf_ready_data.S1_data_mat = ccf_cos_taper_3dim(ccf_ready_data.S1_data_mat,S1_size);
        ccf_ready_data.S2_data_mat = ccf_cos_taper_3dim(ccf_ready_data.S2_data_mat,S2_size);
        disp('Cos-Tapering operation on the signal completed')
   end
    
   if filttype.IsPrefilter
        ccf_ready_data.S1_data_mat = ccf_butterfilt_3dim(ccf_ready_data.S1_data_mat,S1_size,filttype,parameters);
        ccf_ready_data.S2_data_mat = ccf_butterfilt_3dim(ccf_ready_data.S2_data_mat,S2_size,filttype,parameters);
        disp('Butterworth filtering operation on the signal completed')
   end

   if filttype.IsMultiTaper
        [PSI,~] = sleptap(S1_size(3),filttype.NW,filttype.K_taps);
        ccf_ready_data.S1_data_mat = ccf_slepian_multitap_3dim(ccf_ready_data.S1_data_mat,PSI,filttype);
        ccf_ready_data.S2_data_mat = ccf_slepian_multitap_3dim(ccf_ready_data.S2_data_mat,PSI,filttype);
        disp('Multitaper operation on the signal completed'); clear PSI;
   end

   if filttype.IsFTN
        [ b, a ] = get_filter_TFcoeffs(filttype.frange_FTN, parameters.dt);
        ccf_ready_data.S1_data_mat = ccf_FTN_3dim(ccf_ready_data.S1_data_mat,S1_size,b,a);
        ccf_ready_data.S2_data_mat = ccf_FTN_3dim(ccf_ready_data.S2_data_mat,S2_size,b,a);
        disp('FTN operation on the signal completed'); clear b a;
   end

   if filttype.IsOBN
        ccf_ready_data.S1_data_mat = ccf_OBN_3dim(ccf_ready_data.S1_data_mat);
        ccf_ready_data.S2_data_mat = ccf_OBN_3dim(ccf_ready_data.S2_data_mat);
        disp('OBN operation on the signal completed');
   end

   if filttype.IsSpecWhiten
        ccf_ready_data.S1_data_mat = ccf_spectrumwhiten_smooth_3dim(ccf_ready_data.S1_data_mat,S1_size,0.001,filttype);
        ccf_ready_data.S2_data_mat = ccf_spectrumwhiten_smooth_3dim(ccf_ready_data.S2_data_mat,S2_size,0.001,filttype);
        disp('Spectrum Whitening operation on the signal completed');
   end




























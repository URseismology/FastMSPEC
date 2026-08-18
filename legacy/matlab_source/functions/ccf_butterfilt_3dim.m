function [S_data_mat] = ccf_butterfilt_3dim(data,data_dim,filttype,parameters)        
    [bb,aa] = butter(2,filttype.frange_prefilt*2*parameters.dt);
    S_filtfilt = FiltFiltM(bb,aa,reshape(data,[],data_dim(3))')';
    S_data_mat = reshape(S_filtfilt,data_dim(1),data_dim(2),data_dim(3));
    clear S1_filtfilt bb aa
end
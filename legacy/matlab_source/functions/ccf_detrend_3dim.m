function [data_mat] = ccf_detrend_3dim(data,data_dim)
    SZ_detrnd = detrend(reshape(data,[],data_dim(3))')';
    data_mat = reshape(SZ_detrnd,data_dim(1),data_dim(2),data_dim(3));
    clear SZ_detrnd
end

function tapered = ccf_cos_taper_3dim(data,data_dim)
    datalen = data_dim(3);
    M=floor((datalen*5/100)/2+0.5);
    tapered=zeros(1,datalen);
    
    for j=1:datalen
        if j<=M+1
            tapered(j) = 0.5 * ( 1-cos(j*pi/(M+1)));
        elseif (j<datalen - M-1)
            tapered(j) = 1;
        elseif j<=datalen
            tapered(j) = 1 * (0.5 * (1-cos((datalen-j)*pi/(M+1))));
        end
    end
    
    tapered = data .* reshape(tapered, [1, 1, datalen]);
    clear datalen M

addpath('..');
sta1 = 'SA58'; sta2 = 'SA53'; comp='BHZ';
datadir = '../../data/test/raw_data/';
dt=1; Nstart_sec=50; Nstart=Nstart_sec/dt; winlength=4;

list1 = dir([datadir,sta1,'/*',comp,'.sac']);
row_count = 1;
S1_data_mat = [];
S2_data_mat = [];

for ifil = 1:length(list1)
    file1cZ = list1(ifil).name;
    Nchar = length(sta1);
    suffix = file1cZ(Nchar+1:end);
    file2cZ = dir([datadir,sta2,'/',sta2,suffix]);
    if isempty(file2cZ)
        continue
    end

    data1cZ = [datadir,sta1,'/',file1cZ];
    data2cZ = [datadir,sta2,'/',file2cZ.name];

    S1 = readsac(data1cZ);
    S2 = readsac(data2cZ);
    S1Zraw = S1.DATA1; S1Zt = [0:S1.NPTS-1]'*S1.DELTA;
    S2Zraw = S2.DATA1; S2Zt = [0:S2.NPTS-1]'*S2.DELTA;

    if S1.DELTA ~= S2.DELTA
        error('delta mismatch');
    end

    t1num = datenum(S1.NZYEAR,1,S1.NZJDAY,S1.NZHOUR,S1.NZMIN,S1.NZSEC+S1.NZMSEC/1000);
    t2num = datenum(S2.NZYEAR,1,S2.NZJDAY,S2.NZHOUR,S2.NZMIN,S2.NZSEC+S2.NZMSEC/1000);
    offset_sec = (t2num - t1num)*86400;
    if abs(offset_sec) > S1.DELTA
        error('start time mismatch too large');
    end
    S2Zt = S2Zt + offset_sec;

    if (abs(S1.DELTA-dt) >= 0.01*dt) || (abs(S2.DELTA-dt) >= 0.01*dt)
        error('sampling interval mismatch');
    end

    if nnz(S1Zraw) == 0 || nnz(S2Zraw) == 0
        continue
    end
    if length(S1Zraw) < 30000 || length(S2Zraw) < 30000
        continue
    end

    hour_length = winlength;
    nwin = floor(24/hour_length)*2-1;
    win_length = hour_length*60*60/dt;
    last_pt = win_length*.5*(nwin-1)+1+Nstart+win_length;
    if last_pt < length(S1Zraw)
        nwin = nwin + 1;
    end

    for iwin = 1:nwin
        pts_begin = win_length*.5*(iwin-1)+1+Nstart;
        pts_end = pts_begin+win_length;
        if pts_begin > length(S1Zraw) || pts_begin > length(S2Zraw) || pts_end > length(S1Zraw) || pts_end > length(S2Zraw)
            pts_begin = length(S1Zraw)-win_length-Nstart;
            pts_end = pts_begin+win_length;
        end
        tcut = [pts_begin:pts_end] * dt;
        S1Z = interp1(S1Zt, S1Zraw, tcut); S1Z(isnan(S1Z))=0;
        S2Z = interp1(S2Zt, S2Zraw, tcut); S2Z(isnan(S2Z))=0;
        S1_data_mat(row_count,iwin,:) = S1Z;
        S2_data_mat(row_count,iwin,:) = S2Z;
    end
    row_count = row_count + 1;
end

printf('total valid days: %d\n', row_count-1);
printf('S1_data_mat shape: %s\n', mat2str(size(S1_data_mat)));

fftS1Z = fft(S1_data_mat,[],3);
fftS2Z = fft(S2_data_mat,[],3);
coh_trace = fftS2Z .* conj(fftS1Z);
coh_trace = coh_trace ./ abs(fftS1Z) ./ abs(fftS2Z);
coh_trace(isnan(coh_trace)) = 0;
coh_sum = sum(sum(coh_trace,1),2);
coh_sum = reshape(coh_sum,1,[]);
coh_num = size(coh_trace,1)*size(coh_trace,2);
printf('coh_num=%d coh_sum(1:5)=\n', coh_num); disp(coh_sum(1:5))

save('-v7','full_pipeline_octave_out.mat','coh_sum','coh_num');

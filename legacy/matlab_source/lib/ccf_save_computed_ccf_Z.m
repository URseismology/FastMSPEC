function [saved_ccf_path] = ccf_save_computed_ccf_Z(coh_trace,coh_num_size,coh_num,path_vars,parameters,stapairsinfo,stacktype)

    sta1 = stapairsinfo.stanames{1};
    sta2 = stapairsinfo.stanames{2};

    %save full data stack
    if stacktype.IsOutputFullstack
        coh_trace_sum = sum(sum(coh_trace,1),2);
        coh_sum = reshape(coh_trace_sum,1,[]);
        
        ccfZ_fullstack_path = [path_vars.ccf_fullstack_path,'ccf',parameters.strNAMEcomp,'/'];
        save(sprintf('%s%s/%s_%s_f.mat',ccfZ_fullstack_path,sta1,sta1,sta2),'coh_sum','coh_num','stapairsinfo');
        
        saved_ccf_path.fullstack = sprintf('%s%s/%s_%s_f.mat',ccfZ_fullstack_path,sta1,sta1,sta2);
        disp(['Full Stacked Cross Correlation has been computed for ' sta1 '-' sta2  ' sation pairs'])
    end

    %save each day stack
    if stacktype.IsOutputDaystack
        coh_trace_sum = sum(coh_trace,2);
        coh_trace_size = size(coh_trace_sum);
        coh_sum = reshape(coh_trace_sum,coh_trace_size(1),coh_trace_size(3));
        
        ccfZ_daystack_path = [path_vars.ccf_daystack_path,'ccf',parameters.strNAMEcomp,'/'];
        save(sprintf('%s%s/%s_%s_f.mat',ccfZ_daystack_path,sta1,sta1,sta2),'coh_sum','coh_num_size','stapairsinfo');
        
        saved_ccf_path.daystack = sprintf('%s%s/%s_%s_f.mat',ccfZ_daystack_path,sta1,sta1,sta2);
        disp(['Day Stacked Cross Correlation has been computed for ' sta1 '-' sta2  ' sation pairs'])
    end

    %save each month stack
    if stacktype.IsOutputMonthstack
        coh_trace_size = size(coh_trace);
        n_days = coh_trace_size(1); n_data = coh_trace_size(3);
        total_months = ceil(n_days/30);
        coh_trace_month = zeros([total_months,n_data]);
        start_idx=1;end_idx=30;
        month_counter=1;
        for i=1:total_months
            coh_trace_subset = coh_trace(start_idx:end_idx,:,:);
            coh_trace_subset = sum(sum(coh_trace_subset,1),2);
            coh_trace_subset = reshape(coh_trace_subset,1,[]);
            coh_trace_month(month_counter,:) = coh_trace_subset;
            start_idx = end_idx+1;
            end_idx = end_idx + 30;
            if end_idx > n_days
                end_idx = n_days;
            end
        end
        coh_sum = coh_trace_month;
        
        ccfZ_monthstack_path = [path_vars.ccf_monthstack_path,'ccf',parameters.strNAMEcomp,'/'];
        save(sprintf('%s%s/%s_%s_f.mat',ccfZ_monthstack_path,sta1,sta1,sta2),'coh_sum','coh_num_size','stapairsinfo');
        
        saved_ccf_path.monthstack = sprintf('%s%s/%s_%s_f.mat',ccfZ_monthstack_path,sta1,sta1,sta2);
        disp(['Month Stacked Cross Correlation has been computed for ' sta1 '-' sta2  ' sation pairs'])
    end

    clear ccfZ_fullstack_path ccfZ_daystack_path ccfZ_monthstack_path sta1 sta2 coh_sum coh_trace coh_trace_sum
    clear n_days coh_trace_size total_months coh_trace_month start_idx end_idx month_counter coh_trace_subset coh_num_size
    

end
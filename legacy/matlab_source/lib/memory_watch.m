function [totalMB] = memory_watch()
    vars = evalin('caller','whos');
    totalBytes = 0;
    
    for i = 1:length(vars)
        totalBytes = totalBytes + vars(i).bytes;
    end
    
    totalMB = totalBytes / (1024^2);
end
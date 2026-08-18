function tf = contains(str, pattern)
%CONTAINS  Minimal Octave shim for MATLAB's contains(str, pattern),
% single string/pattern case only (as used by mspec_fast.m).
tf = ~isempty(strfind(str, pattern));
end

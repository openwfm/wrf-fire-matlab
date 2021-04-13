function [icl1,icl2]=hzc2icl(hzc,n)
% [icl1,icl2]=hzc2icl(hzc,n)
% In:
%   hzc     vector size 2, coarsening factors
%   n       vector size 3, fine mesh size
% Out:
%   icl1    coarse indices in the horizontal direction 1
%   icl2    coarse indices in the horizontal direction 2
icl1=1:hzc(1):n(1);
if icl1(end) ~= n(1)
    icl1(end+1)=n(1); % extend, end must be coarse
end
icl2=1:hzc(2):n(2);
if icl2(end) ~= n(2)
    icl2(end+1)=n(2);
end
end
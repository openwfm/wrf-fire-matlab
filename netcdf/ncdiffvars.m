function [relerr,ssq,p]=ncdiffvars(s,v1,v2)
% [relerr,ssq,p]=ncdiffvars(s,v1,v2)
% compare two array values v1 and v2 using a description s
% relerr = max relative difference
% ssq = square mean relative difference
% p = p-value, close to +1 -1 means errors significantly in one direction
% for rounding errors expect relerr=const*eps, p small
scale=max(big(v1),big(v2));
relerr=big(v2(:)-v1(:))/(scale+realmin); 
ssq=norm(v2(:)-v1(:))/(max(norm(v1(:)),norm(v2(:)))+realmin); 
d=(v2(:)-v1(:))/(scale+realmin); % scaled diff
avgdiff=mean(d);
stdev=std(d);
n=length(v1(:));
t=sqrt(n)*avgdiff/(stdev+realmin);
p=erf(t);

fprintf('%s max abs %g relative diff %g max %g stdev %g avg %g t-stats %g p-value %g\n',...
    s,scale,relerr,max(d),stdev,avgdiff,t,p) 
end
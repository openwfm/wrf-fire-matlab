function [relerr,ssq,names]=ncdiff(file1,file2,var)
% [relerr,stdev,p]=ncdiff(file1,file2[,'var'])
% compare variable var in 2 netcdf files
% if var is missing, compare all variables
% relerr = max relative difference
% ssq = square mean relative difference
% p = p-value, close to +1 -1 means errors significantly in one direction
% for rounding errors expect relerr=const*eps, p small
if ~exist('var','var'),
    v=ncdump(file1,'-q');
    for i=1:length(v),
        var=v(i).varname;
        try
            [relerr(i),ssq(i),names{i}]=ncdiff(file1,file2,var);
        catch exc
            disp(exc)
            relerr(i)=inf;
        end
    end
    disp('Finished comparing files')
    disp(file1)
    disp(file2)
    return
end
v1=ncread(file1,var);
v2=ncread(file2,var);
s1=size(v1);
s2=size(v2);
if length(s1) ~=  length(s2) | any(s1 ~= s1),
    error([var,' must be same size 1: ',num2str(s1),' 2: ',num2str(s2)]) 
end 
[relerr,ssq]=comp('total',v1,v2);
if relerr>0,
    for i=1:s1(end)    
        s=['slice ',num2str(i)];
        switch length(s1)
            case 4
                comp(s,v1(:,:,:,i),v2(:,:,:,i));
            case 3
                comp(s,v1(:,:,i),v2(:,:,i));
        end
    end
end
names=var;
end

function [relerr,ssq]=comp(s,v1,v2)
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

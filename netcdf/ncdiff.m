function [relerr,ssq,p,name]=ncdiff(file1,file2,var)
% [relerr,ssq,p,name]=ncdiff(file1,file2,var)
% compare variable var in 2 netcdf files
% if var is missing, compare all variables
% return structure s with fiels:
% relerr = max relative difference
% ssq = square mean relative difference
% p = p-value, close to +1 -1 means errors significantly in one direction
% for rounding errors expect relerr=const*eps, p small
% name = variable name
if ~exist('var','var'),
    v=ncdump(file1,'-q');
    for i=1:length(v),
        var=v(i).varname;
        try
            [relerr(i),ssq(i),p(i),name{i}]=ncdiff(file1,file2,var);
        catch exc
             disp(exc)
             relerr(i)=inf;
             p(i)=1;
        end
    end
    disp('Finished comparing files')
    disp(file1)
    disp(file2)
    for i=1:length(p),
       fprintf('%s relative diff %g p-value %g\n',name{i},relerr(i),p(i))
    end
    return
end
v1=ncread(file1,var);
v2=ncread(file2,var);
s1=size(v1);
s2=size(v2);
if length(s1) ~=  length(s2),
    error([var,' must be same numer of dimensions 1: ',num2str(s1),' 2: ',num2str(s2)]) 
end 
if any(s1 ~= s2)
    warning(['using minimum size of 1:', num2str(s1),' 2: ',num2str(s2)])
    s=ones(1,4);
    s(1:length(s1))=min(s1,s2);
    v1=v1(1:s(1),1:s(2),1:s(3),1:s(4));
    v2=v2(1:s(1),1:s(2),1:s(3),1:s(4));
    s1=size(v1);
    s2=size(v2);
end
[relerr,ssq,p]=ncdiffvars('total',v1,v2);
if relerr>0,
    for i=1:s1(end)    
        s=['slice ',num2str(i)];
        switch length(s1)
            case 4
                ncdiffvars(s,v1(:,:,:,i),v2(:,:,:,i));
            case 3
                ncdiffvars(s,v1(:,:,i),v2(:,:,i));
        end
    end
end
name=var;
end

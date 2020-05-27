function var=nclist(filename,q) 
%  info on all variables

quiet=exist('q','var');
fprintf('ncdump: file %s\n',filename);
ncid = netcdf.open(filename,'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdim] = netcdf.inq(ncid);
for varid=1:nvars, % one variable at a time
    var(varid)=ncvarinfo(ncid,varid-1);
    if ~quiet,
        fprintf('%i ',varid);
        dispvarinfo(var(varid));
    end
end
for i=1:nvars % save the native netcdf variable order
    var(i).varid=i-1;
end
[~,ix]=sort({var.varname});
var=var(ix);
netcdf.close(ncid);
end

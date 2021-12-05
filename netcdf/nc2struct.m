function [p,dims]=nc2struct(filename,varnames,gattnames,timestep,p)
% p=nc2struct(filename,varnames,gattnames,timestep,p)
% read from netcdf file to structure
%
% arguments
% input
%   filename        string, name of file to read
%   varnames        cell array of strings, names of variable to read
%   gattnames       cell array of strings, names of global attributes to read
%   times           (optional) matrix of indices in the last dimension to extract (timestep in WRF)
%   p               (optional) the structure to add to
% output
%   p               matlab structure with the specifed variables and attributes as fields
%                   with names in lowercase.  The types are kept and the dimensions 
%                   are not squeezed
%
% example
%   p=nc2struct('wrfinput_d01',{'U','V'},{'DX','DY'})
% will read variables U,V into p.u, p.v and global attributes DX DY into
% p.dx p.dy, respectively

if ~exist('timestep','var'),
    timestep=0;
    t=-1;
    fprintf('nc2struct: reading all timesteps\n')
elseif isscalar(timestep) & isnumeric(timestep),
    fprintf('nc2struct: reading timestep %i only\n',timestep)
    t=timestep-1; % netcdf dimensions start from 0
else
    error('timestep must be numeric scalar')
end
p.timestep=timestep;

fprintf('nc2struct: reading from file %s\n',filename)

try
   ncid = netcdf.open(filename,'NC_NOWRITE');
catch ERR 
   disp(['cannot open NetCDF file ',filename])
   rethrow(ERR);
end
netcdf.close(ncid);

p.filename{1}=filename;

% reading values


for i=1:length({varnames{:}}),
    varname=varnames{i};
    try
        v=ncvar(filename,varname,[]);
    catch ME
        warning(['variable ',varname,' does not exist in file ',filename])
        v=[];
    end
    field = strrep(lower(varname),'-','__'); 
    if ~ isempty(v),
        ndims=length(v.dimlength);
        start=zeros(1,ndims);
        count=v.dimlength;
        if t >= 0, % read only one dimestep
            if(v.dimids(ndims)~=0),
                 warning('id of the last dimension is not 0, is it timestep?')
            end
            start(ndims)=t;
            count(ndims)=1;
        end
        v = ncvar(filename,varname,start, count); 
        p.(field)=double(v.var_value);
        dims.(field)=v.dimlength;
    else 
        p.(field)=[];
    end
end

% reading attributes

for i=1:length({gattnames{:}}),
    gattname=gattnames{i};
    try
        val=ncgetgatt(filename,gattname);
    catch ME
        warning(['global attribute ',gattname,' does not exist in file ',filename])
        val=[];
    end
    p.(lower(gattname))=double(val);
end

end

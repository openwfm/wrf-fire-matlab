function [X,u0]=read_fmw_wrfout(path,timestep)

% if timestep not exist, last timestep
if ~exist('timestep','var'),
    timestep=size(nc2struct(path,{'Times'},{}).times,2);
end
% reading wrfout variables
% p=nc2struct(path,{'U0_FMW','V0_FMW','W0_FMW','HT_FMW','ZSF'},{'DX','DY'},timestep);
p=nc2struct(path,{'U0_FMW','V0_FMW','W0_FMW','HT_FMW','ZSF','U_FMW','W_FMW','V_FMW','UF','VF'},{'DX','DY'},timestep);
% looking at sizes
[nx,ny,nz]=size(p.u0_fmw);
% height profile from the ground
htt_fmw=zeros(nz+1,1);
for z=1:nz
    htt_fmw(z+1)=2*p.ht_fmw(z)-htt_fmw(z);
end
% compute elevation at corners from midpoints
zsf=midpoints2corners(p.zsf);
% creating grid 
[Cx,Cy,Cz]=ndgrid(p.dx*[0:nx],p.dy*[0:ny],htt_fmw);
% adding ZSF
for z=1:nz+1
    Cz(:,:,z)=Cz(:,:,z)+zsf;
end
% generate and check mesh
X={Cx,Cy,Cz};
check_mesh(X);
% compute initial wind at the corners from midpoints
%Xu0=midpoints2corners(p.u0_fmw);
%Yu0=midpoints2corners(p.v0_fmw);
%Zu0=midpoints2corners(p.w0_fmw);
% generate initial wind at corners of the cells
%u0={Xu0,Yu0,Zu0};
u0={p.u0_fmw,p.v0_fmw,p.w0_fmw};

end

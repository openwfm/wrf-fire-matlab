function [X,u0]=read_fmw_wrfout(path,timestep)
% if not specified, first time step
if ~exist('timestep','var'),
    timestep=1
% reading wrfout variables
p=nc2struct(path,{'U0_FMW','V0_FMW','W0_FMW','HT_FMW','ZSF'},{'DX','DY'},timestep)
% looking at horizontal size
[nx,ny,nz]=size(p.u0_fmw);
% elevation at corners
ht_fmw=[0;(p.ht_fmw(1:end-1)+p.ht_fmw(2:end))/2];
ht_fmw=[ht_fmw;ht_fmw(end)+1.2810*(ht_fmw(end)-ht_fmw(end-1))];
zsf=zeros(nx+1,ny+1);
zsf(2:end-1,2:end-1)=(p.zsf(1:end-1,1:end-1)+p.zsf(1:end-1,2:end)+p.zsf(2:end,1:end-1)+p.zsf(2:end,2:end))/4;
% creating grid 
[Cx,Cy,Cz]=ndgrid(p.dx*[0:nx],p.dy*[0:ny],ht_fmw);
% adding ZSF
for z=1:nz+1
    Cz(:,:,z)=Cz(:,:,z)+zsf;
end
% compute wind at the corners from midpoints
%p.u0_fmw
%p.v0_fmw
%p.w0_fmw
% generate mesh and wind cells
X={Cx,Cy,Cz};
check_mesh(X);
u0={Xu0,Yu0,Zu0};
end

function p=read_fmw_3d(path,frame)
% p=read_fmw_3d(path,frame)
% read initial and resulting wind at cell centers and compute their coordinates
p=nc2struct(path,{'U0_FMW','V0_FMW','W0_FMW','HT_FMW','ZSF','FWH',...
    'U_FMW','W_FMW','V_FMW','UF','VF','FWH','U','V','W','ZS','HGT',...
    'PHB','PH','P','PB'},{'DX','DY'},frame);
if ~isempty(p.u0_fmw)
    disp('fmw variables, computing derived')
else
    return
end
% WRF atm wind at cell centers (theta points) 
% see https://wiki.openwfm.org/wiki/How_to_interpret_WRF_variables
% cell centers coordinates
[nx,ny,nz]=size(p.w);
nz=nz-1;
p.xc=zeros(nx,ny,nz);  
p.yc=zeros(nx,ny,nz);
p.zc=zeros(nx,ny,nz);
ph = p.phb + p.ph;
p.zw=ph/9.81;
p.zc =0.5*(ph(:,:,1:end-1)+ph(:,:,2:end))/9.81
[p.xc(:,:,1),p.yc(:,:,1)]=ndgrid(p.dx*([1:nx]-0.5),p.dy*([1:ny]-0.5));
for k=2:nz
    p.xc(:,:,k) = p.xc(:,:,1);
    p.yc(:,:,k) = p.yc(:,:,1);
end
% average the wind from staggered grid to cell centers
p.uc = 0.5*(p.u(1:end-1,:,:) + p.u(2:end,:,:));
p.vc = 0.5*(p.v(:,1:end-1,:) + p.v(:,2:end,:));
p.wc = 0.5*(p.w(:,:,1:end-1) + p.w(:,:,2:end));
p.uvwc=sqrt(p.uc.^2+p.vc.^2+p.wc.^2); % 3D velocity
p.uvc=sqrt(p.uc.^2+p.vc.^2); % horizontal velocity
%*** fmw arrays
% fmw cell center coordinates, same as wind
[nxf,nyf,nzf]=size(p.u0_fmw);
% fire step sizes
p.dxf=p.dx/(nxf/nx);
p.dyf=p.dx/(nxf/nx);
% coordinates of centers
p.xcf=zeros(nxf,nyf,nzf);  
p.ycf=zeros(nxf,nyf,nzf);
p.zcf=zeros(nxf,nyf,nzf);
[p.xcf(:,:,1),p.ycf(:,:,1)]=ndgrid(p.dxf*([1:nxf]-0.5),p.dyf*([1:nyf]-0.5));
ht0=[0;p.ht_fmw]; p.htcf=(ht0(2:end)+ht0(1:end-1))/2;
p.zcf(:,:,1) = p.zsf + p.htcf(1);
for k=2:nzf
    p.xcf(:,:,k) = p.xcf(:,:,1);
    p.ycf(:,:,k) = p.ycf(:,:,1);
    p.zcf(:,:,k) = p.zsf + p.htcf(k);
end
%***
% fmw wind velocities
p.uv0=sqrt(p.u0_fmw.^2+p.v0_fmw.^2);
p.uvw0=sqrt(p.u0_fmw.^2+p.v0_fmw.^2+p.w0_fmw.^2);
p.uvf=sqrt(p.u_fmw.^2+p.v_fmw.^2);
p.uvwf=sqrt(p.u_fmw.^2+p.v_fmw.^2+p.w_fmw.^2);
% interpolated to fire wind height fwh
p.ufvf=sqrt(p.uf.^2+p.vf.^2);
end

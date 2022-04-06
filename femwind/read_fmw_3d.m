function p=read_fmw_3d(path)
% p=read_fmw_3d(path)
% read intial and resulting wind at cell centers and compute their coordinates
p=nc2struct(path,{'U0_FMW','V0_FMW','W0_FMW','HT_FMW','ZSF',...
    'U_FMW','W_FMW','V_FMW','UF','VF','FWH','U','V','W'},{'DX','DY'},1);
[nxf,nyf,nzf]=size(p.u0_fmw);
p.xcf=zeros(nxf,nyf,nzf);  % cell center coordinates, same as wind
p.ycf=zeros(nxf,nyf,nzf);
p.zcf=zeros(nxf,nyf,nzf);
[p.xcf(:,:,1),p.ycf(:,:,1)]=ndgrid(p.dx*[1:nxf],p.dy*[1:nyf]);
ht0=[0;p.ht_fmw]; p.htcf=(ht0(2:end)+ht0(1:end-1))/2;
p.zcf(:,:,1) = p.zsf + p.htcf(1);
for k=2:nzf
    p.xcf(:,:,k) = p.xcf(:,:,1);
    p.ycf(:,:,k) = p.ycf(:,:,1);
    p.zcf(:,:,k) = p.zsf + p.htcf(k);
end
p.uv0=sqrt(p.u0_fmw.^2+p.v0_fmw.^2);
p.uvw0=sqrt(p.u0_fmw.^2+p.v0_fmw.^2+p.w0_fmw.^2);
p.uv=sqrt(p.u_fmw.^2+p.v_fmw.^2);
p.uvw=sqrt(p.u_fmw.^2+p.v_fmw.^2+p.w_fmw.^2);
p.uvf=sqrt(p.uf.^2+p.vf.^2);
end

function p=read_fmw_3d(path)
% p=read_fmw_3d(path)
% read intial and resulting wind at cell centers and compute their coordinates
p=nc2struct(path,{'U0_FMW','V0_FMW','W0_FMW','HT_FMW','ZSF',...
    'U_FMW','W_FMW','V_FMW','UF','VF','FWH'},{'DX','DY'},1);
[nx,ny,nz]=size(p.u0_fmw);
p.xc=zeros(nx,ny,nz);  % cell center coordinates, same as wind
p.yc=zeros(nx,ny,nz);
p.zc=zeros(nx,ny,nz);
[p.xc(:,:,1),p.yc(:,:,1)]=ndgrid(p.dx*[1:nx],p.dy*[1:ny]);
p.ht0=[0;p.ht_fmw]; p.htc=(p.ht0(2:end)+p.ht0(1:end-1))/2;
p.zc(:,:,1) = p.zsf + p.htc(1);
for k=2:nz
    p.xc(:,:,k) = p.xc(:,:,1);
    p.yc(:,:,k) = p.yc(:,:,1);
    p.zc(:,:,k) = p.zsf + p.htc(k);
end
p.uv0=sqrt(p.u0_fmw.^2+p.v0_fmw.^2);
p.uvw0=sqrt(p.u0_fmw.^2+p.v0_fmw.^2+p.w0_fmw.^2);
p.uv=sqrt(p.u_fmw.^2+p.v_fmw.^2);
p.uvw=sqrt(p.u_fmw.^2+p.v_fmw.^2+p.w_fmw.^2);
p.uvf=sqrt(p.uf.^2+p.vf.^2);
end

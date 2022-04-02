function plot_fmw_3d(path,level,scale,stride)
p=nc2struct(path,{'U0_FMW','V0_FMW','W0_FMW','HT_FMW','ZSF',...
    'U_FMW','W_FMW','V_FMW','UF','VF','FWH'},{'DX','DY'},1);
[nx,ny,nz]=size(p.u0_fmw);
xc=zeros(nx,ny,nz);  % cell center coordinates, same as wind
yc=zeros(nx,ny,nz);
zc=zeros(nx,ny,nz);
[xc(:,:,1),yc(:,:,1)]=ndgrid(p.dx*[1:nx],p.dy*[1:ny]);
ht=[0;p.ht_fmw]; htc=(ht(2:end)+ht(1:end-1))/2;
zc(:,:,1) = p.zsf + ht(1);
for k=2:nz
    xc(:,:,k) = xc(:,:,1);
    yc(:,:,k) = yc(:,:,1);
    zc(:,:,k) = p.zsf + ht(k);
end
figure(1)
plot_wind_3d({xc,yc,zc},{p.u0_fmw,p.v0_fmw,p.w0_fmw},level,scale,stride)
title('Initial wind')
xlabel('m');ylabel('m');zlabel('m')
figure(2)
plot_wind_3d({xc,yc,zc},{p.u_fmw,p.v_fmw,p.w_fmw},level,scale,stride)
title('Mass consistent wind')
xlabel('m');ylabel('m');zlabel('m')


    function vprof(i,j)
        % add vertical wind profile at location i,j
        ,
    end
end

function p=plot_fmw_3d(path,level,scale,stride)
p=read_fmw_3d(path)
figure(1)
plot_wind_3d({p.xcf,p.ycf,p.zcf},{p.u0_fmw,p.v0_fmw,p.w0_fmw},level,scale,stride)
title('Initial wind')
xlabel('m');ylabel('m');zlabel('m')
figure(2)
plot_wind_3d({p.xcf,p.ycf,p.zcf},{p.u_fmw,p.v_fmw,p.w_fmw},level,scale,stride)
title('Mass consistent wind')
xlabel('m');ylabel('m');zlabel('m')
    function vprof(i,j)
        % add vertical wind profile at location i,
    end
end

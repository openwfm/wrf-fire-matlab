function err=prol_rest_err(hzc,icl3,X,params)
% err=prol_rest_err(hzc,icl3,X,params)
% test prolongation and restriction
% err=0,return
n = size(X{1});
[icl1,icl2]=hzc2icl(hzc,n);
nc = [length(icl1),length(icl2),length(icl3)];
P =coarsening_P(hzc,icl3,X,params);
uc = rand(nc);
u0 = reshape(P*uc(:),n);
u1 = prolongation(uc,hzc,icl3,X,params);
err1 = big(u0-u1);
tol = eps(single(1))*10*big(u0);
if err1 > tol 
    err1,tol
    warning('prol_restr_err: prolongation error vs. matrix P too large')
end
if exist('fortran/prolongation_test.exe')
    disp('testing if same result in fortran')
    u2 = prolongation_fortran(uc,hzc,icl3,X,params);
    err2 = big(u0-u2);
    if err2 > tol
    err2,tol
        warning('prol_restr_err: prolongation error fortran too large')
    end
end

u = rand(n);
uc0 = reshape(P'*u(:),nc);
uc1 = restriction(u,hzc,icl3,X,params);
errc = big(uc0-uc1);
if errc > tol
    errc,tol
    warning('prol_restr_err: restriction error vs. matrix P too large')
end
err=max([err1,err2,errc]);
if err>tol,
    error('prol_rest_err: error too large')
else
    fprintf('prol_rest_err %g OK\n',err)
end

end

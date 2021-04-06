function err=prol_rest_err(hzc,icl3,X,params)
% err=prolongation_err(icl,X,params)
n = size(X{1});
[icl1,icl2]=hzc2icl(hzc,n);
nc = [length(icl1),length(icl2),length(icl3)];
tol = eps*10   
P =coarsening_P(hzc,icl3,X,params);
uc = rand(nc);
u0 = reshape(P*uc(:),n);
u1 = prolongation(uc,hzc,icl3,X,params);
err = big(u0-u1)
if err > tol*big(u0)
    error('prol_restr_err: prolongation error too large')
end
u = rand(n);
uc0 = reshape(P'*u(:),nc);
uc1 = restriction(u,hzc,icl3,X,params);
errc = big(uc0-uc1)
if errc > tol*big(uc0)
    error('prol_restr_err: restriction error too large')
end
end

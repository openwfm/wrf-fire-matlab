function err=prol_rest_err(icl,X,params)
% err=prolongation_err(icl,X,params)
n = size(X{1});
nc = cellfun(@length,icl);
P =coarsening_P(icl,X,params);
uc = rand(nc);
u0 = reshape(P*uc(:),n);
u1 = prolongation(uc,icl,X,params);
err = big(u0-u1)
if err > 1e-14*big(u0)
    error('prol_restr_err: prolongation error too large')
end
u = rand(n);
uc0 = reshape(P'*u(:),nc);
uc1 = restriction(u,icl,X,params);
errc = big(uc0-uc1)
if errc > 1e-14*big(uc0)
    error('prol_restr_err: restriction error too large')
end
end

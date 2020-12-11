function condition=scond(K)
R= chol(K);
afun = @(v) R\(R'\v)
emax= eigs(K,1,'largestreal')
emin= 1/eigs(afun,size(K,1),1,'largestreal','IsFunctionSymmetric',true)
if size(K,1)<5000,
    e = eig(full(K)); emax_direct=max(e), emin_direct=min(e)
end
condition= emax/emin;
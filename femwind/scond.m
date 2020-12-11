function condition=scond(K)
relative_nonsymmetry=big(K-K')/big(K)
Ks = (K+K')/2;
e= eigs(Ks,2,'bothendsreal')
condition= max(e)/min(e);

function P=cell_P(A,B,D)

X = B'*A*B;
L = chol(X,'lower');
U = L\eye(size(X));
ib = L'\U;
X = D*ib*D';
L = chol(X,'lower');
U = L\eye(size(X));
ibd = L'\U;
P=ib - ib*D'*ibd*D*ib;

end
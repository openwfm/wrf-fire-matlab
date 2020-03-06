function P=cell_P(A,B,D)

X = inv(B'*A*B);
P = X - X*D'*((D*X*D')\(D*X));

end
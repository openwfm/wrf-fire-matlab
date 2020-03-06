function cell_P(A,B,D)

Bt=B.';
Dt=D.';
ib=inv(Bt*A*B);
ibd=inv(D*ib*Dt);
P=ib - ib*Dt*ibd*D*ib;

end
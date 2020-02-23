% matrix implementation of the saddlepoint method

% repeated subexpressions
Bt=B.';
Dt=D.';
Ct=C.';
ib=inv(Bt*A*B);
ibd=inv(D*ib*Dt);

% solution matrices
P=ib - ib*Dt*ibd*D*ib;
M =C*P*Ct;
rhs=C*P*Bt*A*v0;

% solution components
q=M\rhs;
v=P*(Bt*A*v0-Ct*q);
p=ibd*D*ib*(Bt*A*v0-Ct*q);

% test
res=[Bt*A*B*v + Dt*p + Ct*q - Bt*A*v0
    D*v
    C*v];
err=norm(res,inf)
disp('see matrix M of the system solved')

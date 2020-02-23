disp('run a setup first!')
% repeated subexpressions
Bt=B.';
Dt=D.';
Ct=C.';
ib=inv(Bt*B);
ibd=inv(D*ib*Dt);

% solution matrices
P=ib - ib*Dt*ibd*D*ib;
M =C*P*Ct;
rhs=(C*P*Bt)*v0;

% solution components
q=M\rhs;
v=P*(Bt*v0-Ct*q);
p=ibd*D*ib*(Bt*v0-Ct*q);

% test
res=[Bt*B*v + Dt*p + Ct*q - Bt*v0
    D*v
    C*v];
err=norm(res,inf)
disp('see matrix M for the system to be solved')

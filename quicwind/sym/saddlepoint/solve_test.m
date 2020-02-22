% repeated subexpressions
setup
Bt=B.';
Dt=D.';
Ct=C.';
ib=inv(Bt*B)
ibd=inv(D*ib*Dt)
P=ib - ib*Dt*ibd*D*ib
% solution
M =C*P*Ct
rhs=(C*P*Bt)*v0
q=M\rhs
v=P*(Bt*v0-Ct*q)
% aux 
p=ibd*D*ib*(Bt*v0-Ct*q)
% verify
res = [Bt*B*v + Dt*p + Ct*q - v0
       D*v
       C*v
       ]
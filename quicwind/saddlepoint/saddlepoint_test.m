B=zeros(9);
for i=1:3:7
    B(i:i+2,i:i+2)=rand(3);
end
D = [1 1 1 0 0 0 0 0 0 
     0 0 0 1 1 1 0 0 0
     0 0 0 0 0 0 1 1 1];
% D = rand(2,4
C = [1 0 0 1 0 0 0 0 0
     0 0 0 1 0 0 0 1 0    
     0 0 0 0 0 0 1 0 0]
v0 = rand(9,1)

% repeated subexpressions
Bt=B.';
Dt=D.';
Ct=C.';
ib=inv(Bt*B)
ibd=inv(D*ib*Dt)

% solution matrices
P=ib - ib*Dt*ibd*D*ib
M =C*P*Ct
rhs=(C*P*Bt)*v0

% solution components
q=M\rhs
v=P*(Bt*v0-Ct*q)
p=ibd*D*ib*(Bt*v0-Ct*q)

% test
res=[Bt*B*v + Dt*p + Ct*q - Bt*v0
    D*v
    C*v]
err=norm(res,inf)

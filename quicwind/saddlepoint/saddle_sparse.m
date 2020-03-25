% matrix implementation of the saddlepoint method

% repeated subexpressions
Bt=B.';
Dt=D.';
Ct=C.';
t=tic
% ib_o=inv(Bt*A*B);
% ibd_o=inv(D*ib_o*Dt);
disp('computing inverses')
U=chol(Bt*A*B);
ib=U\(U'\speye(size(U)));
U=chol(D*ib*Dt);
ibd=U\(U'\speye(size(U)));
% err=big(ibd-ibd_o)
disp(num2str(toc(t))),t=tic;

% solution matrices
disp('setting up equations')
P=ib - ib*Dt*ibd*D*ib;
M =C*P*Ct;
rhs=C*P*Bt*A*v0;
disp(num2str(toc(t))),t=tic;

disp('solving the equations')
% solution components
q=M\rhs;
v=P*(Bt*A*v0-Ct*q);
p=ibd*D*ib*(Bt*A*v0-Ct*q);
disp(num2str(toc(t))),t=tic;

% test
res=[Bt*A*B*v + Dt*p + Ct*q - Bt*A*v0
    D*v
    C*v];
err=norm(res,inf)

function g=gradbfs(x)  
% gradients of basis functions
% in:
%   x   vector size 3
% out:
%   g   matrix with each row gradient of one basis function
    Nb=8;
    ib =[  % coordinates of basis functions
    -1    -1    -1
     1    -1    -1
    -1     1    -1
     1     1    -1
    -1    -1     1
     1    -1     1
    -1     1     1
     1     1     1];
    g=zeros(Nb,3);
    for k=1:Nb
        % b(k) = (1+ib(k,1)*x(1))*(1+ib(k,2)*x(2))*(1+ib(k,3)*x(3)); % basis functions
        g(k,:)= [ib(k,1)*(1+ib(k,2)*x(2))*(1+ib(k,3)*x(3)),... % d/d(x(1))
                (1+ib(k,1)*x(1))*ib(k,2)*(1+ib(k,3)*x(3)),... % d/d(x(2))
                (1+ib(k,1)*x(1))*(1+ib(k,2)*x(2))*ib(k,3)]/8; % d/d(x(3))
    end
end
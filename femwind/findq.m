function q=findq(a,n)
% q=findq(a,n,tol)
% find q such that 1+q+q^2 +...+q^(n-1)=a if a>1 else return q=1
% in:
%   a,n size 1
% out:
%   q size 1
if a>n,
    q=1.01;
    for k=1:10000,
        qnew = (1+a*(q-1))^(1/n);
        diff = qnew - q;
        q = qnew;
        % fprintf('iteration %i q=%g diff=%g\n',k,q,diff)
        if abs(diff)<eps
            break
        end
    end
else
    q=1;
end

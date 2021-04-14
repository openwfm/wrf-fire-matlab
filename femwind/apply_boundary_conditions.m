function [Kb,Fb]=apply_boundary_conditions(K,F,X)
% apply zero dirichet boundary condition on sides and on top
% in:
%   K   stiffness matrix
%   F   load vector
%   X   X{1} has the size of the domain
% out:
%   Kb   updated
%   Fb   updated

n = size(X{1});
nn=prod(n);
[i1,i2,i3]=ind2sub(n,1:nn);
bc = i1==1 | i1 ==n(1) | i2 ==1 | i2 == n(2) | i3 == n(3);  % sides and top only, bottom left alone
bx = sub2ind(n,i1(bc),i2(bc),i3(bc)); % matrix indices with zero boundary condition
% check number of nodes with boundary condition
nbc = (n(1)-2)*(n(2)-2)+2*(n(1)-1)*n(3)+2*(n(2)-1)*n(3);
if nnz(bc) ~= nbc,
    error('wrong number of boundary nodes')
end
Kb=[]; Fb=[];
if ~isempty(K)
    if any(size(K)~=nn)
        error('apply_boundary_conditions: inconsistent size K')  
    end
    Kb=K;
    Kb(bx,:)=0;
    Kb(:,bx)=0;
    Kb(bx,:)=0;
    Kb(bx,bx)=big(diag(K))*speye(length(bx));
end
if ~isempty(F)
    Fb=F;
    Fb(bx)=0;
end
end

function check_nonzeros(K,X)
% check_nonzeros(K,X)
% test if matrix K has structure consistent with finite element hexa grid 
% with coordinates X
n=size(X{1});
nn=prod(n);
if any(size(K)~=nn)
    error('K size inconsistent with mesh X')
end
if any(any((K==0)~=(K'==0)))
    error('K does not have symmetric structure')
end
for j=1:nn
    [j1,j2,j3]=ind2sub(n,j);
    ii=find(K(:,j));
    [i1,i2,i3]=ind2sub(n,ii);
    d=max([abs(i1-j1),abs(i2-j2),abs(i3-j3)],[],2);
    f=find(d>1);
    if f
        for xj=f(:)'
            j=jj(xj);
            fprintf('%i mesh node %i %i %i and %i mesh node %i %i %i distance %i value %g\n',...
                i,i1g,i2,i3,j,j1(xj),j2(xj),j3(xj),d(xj),K(i,j))
        end
    end
end
end
            

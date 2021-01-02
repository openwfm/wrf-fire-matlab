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
    [ii,jj,aij]=find(K);
    [i1,i2,i3]=ind2sub(n,ii);
    [j1,j2,j3]=ind2sub(n,jj);
    d=max([abs(i1-j1),abs(i2-j2),abs(i3-j3)],[],2);
    f=find(d>1);
    if f
        for xj=f(:)'
            fprintf('%i mesh node %i %i %i and %i mesh node %i %i %i distance %i value %g\n',...
                ii(xj),i1(xj),i2(xj),i3(xj),...
                jj(xj),j1(xj),j2(xj),j3(xj),d(xj),aij(xj))
        end
        error('K inconsistent with hexa structure nonzeros at distance more than 1')
    end
end
            

function check_nonzeros(Kc,Xc,P,K,X)
% check_nonzeros(Kc,Xc,P,K,X)
% checking nonzeros structure consistent with hexa grid
    fprintf('matrix size %g nonzeros %g density %g%%\n',...
        length(Kc),nnz(Kc),100*nnz(Kc)/length(Kc)^2)
    disp('check_nonzeros: check if structure consistent with hexa grid')
    nc=size(Xc{1});
    nnc=prod(nc);
    if any(size(Kc)~=nnc)
        error('Kc size inconsistent with mesh Xc')
    end
    if any(any((Kc==0)~=(Kc'==0)))
        error('K does not have symmetric structure')
    end
    [ii,jj,aij]=find(Kc);
    [i1,i2,i3]=ind2sub(nc,ii);
    [j1,j2,j3]=ind2sub(nc,jj);
    d=max([abs(i1-j1),abs(i2-j2),abs(i3-j3)],[],2);
    f=find(d>1);
    if f
        for xj=f(:)'
            i = ii(xj);
            j = jj(xj);
            fprintf('%i mesh node %i %i %i and %i mesh node %i %i %i distance %i value %g\n',...
                i,i1(xj),i2(xj),i3(xj),...
                j,j1(xj),j2(xj),j3(xj),d(xj),aij(xj))
            if exist('P','var') & exist('X','var')
                n = size(X{1});
                fi = find(P(:,i))';
                [fi1,fi2,fi3]=ind2sub(n,fi);
                [i1,i2,i3]=ind2sub(n,i);
                fj = find(P(:,j))';
                [fj1,fj2,fj3]=ind2sub(n,fj);
                [j1,j2,j3]=ind2sub(n,j);
                fprintf('coarse %i at %i %i %i interpolates to\n',i,i1,i2,i3)
                for k=1:length(fi)
                    fprintf('fine   %i at %i %i %i\n',fi(k),fi1(k),fi2(k),fi3(k))
                end
                fprintf('coarse %i at %i %i %i interpolates to\n',j,j1,j2,j3)
                for k=1:length(fj)
                    fprintf('fine   %i at %i %i %i\n',fj(k),fj1(k),fj2(k),fj3(k))
                end
                disp('connected by nonzeros of fine matrix')
                [fii,fjj,fval]=find(K(fi,fj));
                sub=sparse(fi(fii),fj(fjj),fval);
                disp(sub)
                error('K inconsistent with hexa structure nonzeros at distance more than 1')
            end
        end
    end
end
            

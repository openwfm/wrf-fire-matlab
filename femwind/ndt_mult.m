function y=ndt_mult(K,x,doprint)
%y=nd_mult(K,x,doprint)
% multiply vector x by matrix K from nd_assembly
if ~exist('doprint','var')
    doprint=0;
end
[n1,n2,n3,m]=size(K);
t = ndt_storage_table(m); 
u=reshape(x,n1,n2,n3);  % make x into grid vector if needed
y=zeros(n1,n2,n3);
% loop in WRF ordering (j,k,i)
for i2=1:n2 
    for i3=1:n3   
        for i1=1:n1
            % global index of row = node i
            for j3=-1:1   % relative index j of neighbor = nonzero in the row i               
                for j2=-1:1               
                    for j1=-1:1
                        % global index triple of neighbor node i+j out of bounds? 
                        if ~(i1+j1<1 || i1+j1>n1 || i2+j2<1 || i2+j2>n2 || i3+j3<1 || i3+j3>n3 )
                            % offset m of the m where the entry (i,i+j) is stored, 0 if this row 
                            % in fortran we won't have the 2+ because
                            % the array t will be indexed -1:1
                            m3 = t(3,2+j1,2+j2,2+j3);
                            m2 = t(2,2+j1,2+j2,2+j3);
                            m1 = t(1,2+j1,2+j2,2+j3);
                            % index of the matrix entry (i,i+j) in K(i+m,:) 
                            jx = t(4,2+j1,2+j2,2+j3);
                            % contribution of K(i,i+j)*x(i+j) if index not out of bounds
                            y(i1,i2,i3)=y(i1,i2,i3)+K(i1+m1,i2+m2,i3+m3,jx)*u(i1+j1,i2+j2,i3+j3);
                            switch doprint
                                case 1   % matlab ordering
                                    fprintf('i=(%g %g %g) j=(%g %g %g) m=(%g %g %g) jx=%g i+m=(%g %g %g) i+j=(%g %g %g)\n',...
                                           [ i1,i2,i3,    j1,j2,j3,    m1,m2,m3,  jx, i1+m1,i2+m2,i3+m3, i1+j1,i2+j2,i3+j3 ])
                                case 2   % WRF storage ordering
                                    fprintf('ikj: row=(%g %g %g) rel.column=(%g %g %g) rel.row.stored=(%g %g %g) at %g\n',...
                                                    [ i1,i3,i2,             j1,j3,j2,             m1,m3,m2,    jx])
                                case 3   % print code
                                    if i1 == 2 && i2 == 2 && i3 == 2 
                                        fprintf('        kmat(i%s,k%s,j%s,%2i)*u(i%s,k%s,j%s) +  &\n',...
                                            pm(m1),pm(m3),pm(m2),jx,pm(j1),pm(j3),pm(j2))
                                    end
                            end 
                        end
                    end
                 end
            end
        end
    end
end
end

function s=pm(j)
switch j
    case 1 
        s='+1';
    case -1
        s='-1';
    case 0
        s='  ';
    otherwise
        error('pm: argument must be in -1:1')
end
end

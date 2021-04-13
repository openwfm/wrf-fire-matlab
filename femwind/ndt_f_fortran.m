function err=ndt_f_fortran(A,X,u0,iflags)

% compare the result of matlab and fortran
y= ndt_f_assembly(A,X,u0);

write_array_nd(swap23(X{1}),'X');
write_array_nd(swap23(X{2}),'Y');
write_array_nd(swap23(X{3}),'Z');

write_array(A,'A')
write_array(iflags,'iflags')

system('./fortran/ndt_f_test.exe');

y_f=read_array_nd('Fvec');

err= norm(y(:)-y_f(:),inf);


function u=swap23(v)
[n1,n2,n3,n4]=size(v);
u=zeros(n1,n3,n2,n4);
for l=1:n4
    for k=1:n3
        for j=1:n2
            for i=1:n1
                 u(i,k,j,l)=v(i,j,k,l);
            end
        end
    end
end
end


end

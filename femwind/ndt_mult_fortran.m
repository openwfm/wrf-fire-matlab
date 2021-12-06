function err=ndt_mult_fortran(kmat,u)

y   = ndt_mult(kmat,u);

% test vs sparse first
kmat_s = ndt_convert(kmat,'sparse');
y_s = kmat_s*u(:);
err = big(y_s - y(:));
if err > eps*10*big(u)*big(kmat)
    err, error('error ndt multiply vs sparse too large')
end
    
% test vs fortran ndt_mult

write_array_nd(swap23(kmat),'kmat');
write_array_nd(swap23(u),'u');

system('./fortran/ndt_mult_test.exe');

y_f=swap23(read_array_nd('y'));

err=norm(y(:)-y_f(:),inf);

end

function Kb=ndt_boundary_conditions_fortran(K,params); 
% Kb=ndt_boundary_conditions_fortran(K,params); 
% call fortran version and compare results
% in:
%     K global stiffness matrix ndt14 format

[n1,n2,n3,st] = size(K);
if st~=14, 
    error('K must be in ndt 14 storage scheme')
end

%Writing to file for use by fortran tester
write_array_nd(swap23(K),'K');

system('./fortran/ndt_boundary_conditions_test.exe');

Kb_f = swap23(read_array_nd('Kb')); % boundary condition from fortran

% the same in matlab
Ks = ndt_convert(K,'sparse');
[Kb_ms,~]=apply_boundary_conditions(Ks,[],{K(:,:,:,1)});

% convert matlab output to nd14 and compare
Kb_m =  ndt_convert(Kb_ms,n1,n2,n3,st);
err = big(Kb_m - Kb_f)

% convert fortran output to sparse and compare
Kb_fs = ndt_convert(Kb_f,'sparse');
err_s = big(Kb_ms - Kb_fs)
if err_s > eps(single(1.0))*big(Kb_fs)*10
    error('error too large')
end
end


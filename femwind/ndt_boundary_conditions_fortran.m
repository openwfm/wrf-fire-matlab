function Kb=ndt_boundary_conditions_fortran(K,params); 
% Kb=ndt_boundary_conditions_fortran(K,params); 
% call fortran version and compare results
% in:
%     K global stiffness matrix ndt14 format

m=size(K,4);
if m~=14, 
    error('K must be in ndt 14 storage scheme')
end

%Writing to file for use by fortran tester
write_array_nd(swap23(K),'K');

system('./fortran/ndt_boundary_conditions_test.exe');

Kb = swap23(read_array_nd('Kb'));

% the same in matlab
Ks = ndt_convert(K,'sparse');
[Ksb,~]=apply_boundary_conditions(Ks,[],{K(:,:,:,1)});

% convert to sparse and compare
Kbs = ndt_convert(Kb,'sparse');
err = big(Ksb - Kbs)
if err > eps(single(1.0))*big(Ksb)*10
    error('error too large')
end
end


function Fb=vec_boundary_conditions_fortran(F,params); 
% Fb=vec_boundary_conditions_fortran(F,params); 
% call fortran version and compare results
% in:
%     F mesh vector 

%Writing to file for use by fortran tester
write_array(swap23(F),'F');

system('./fortran/vec_boundary_conditions_test.exe');

Fb = swap23(read_array('Fb'));

% the same in matlab
[~,Fbm]=apply_boundary_conditions([],F,{F});

err = big(Fbm - Fb)
if err > eps(single(1.0))*big(Fbm)*10
    error('error too large')
end
end


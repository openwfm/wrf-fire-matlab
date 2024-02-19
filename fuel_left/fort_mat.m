function fort_mat(s,x)
%s is string with variable name
%x is matrix to be converted into something that can be cut/paste into
%fortran code

[n,m] = size(x);
format long
%loop through and print
for i = 1:n
    for j = 1:m
        fprintf('%s(%d,%d) = %3.16f \n',s,i,j,x(i,j));
    end
end
fprintf('\n')



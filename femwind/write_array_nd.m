function write_array_nd(a,name)
% write_array_nd(a,name)
% Purpose: write nd matrix 
%
% Arguments
% a         nd matrix, the array to be written
% name  string, the name of the file
% 
if ~isnumeric(a)
    error('matrix must be numeric array')
end
filename=[name,'.txt'];
s = size(a);
n = length(s);
fprintf('writing matrix size %g %g %g %g %g %g %g',s) 
fprintf('to file %s\n',filename)
h=fopen(filename,'w');
fprintf(h,'%i\n',456,n,s);
fprintf(h,'%.12g\n',a);
fclose(h);
end

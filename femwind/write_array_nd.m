function write_array_nd(a,filename)
% write_array_nd(a,filename)
% Purpose: write nd matrix 
%
% Arguments
% filename  string, the name of the file
% a         nd matrix, the array to be written
% 
if ~isnumeric(a)
    error('matrix must be numeric array')
end
s = size(a);
n = length(s);
h=fopen(filename,'w');
fprintf(h,'%i\n',981,n,s);
fprintf(h,'%.12g\n',a);
fclose(h);
end

function a=read_array_nd(name)
% a=read_arrray_nd(filename)
% read nd matrix written by write_array_nd
filename = [name,'.txt'];
fprintf(['reading matrix from file ',filename])
d=load(filename);
if d(1) ~= 456,
    error('not written by write_array_nd')
end
n=d(2);     % number of dimensions
bs=3;       % beginning of size
es=bs+n-1;  % index to to end of size
s=d(bs:es); % array size
fprintf(' size %g %g %g %g %g %g %g',s)
fprintf('\n')
if length(d)~=prod(s) + es
    error('wrong number of terms in the file')
end
a=d(es+1:end);
if length(a)>1
    a=reshape(d(es+1:end),s');
end
end


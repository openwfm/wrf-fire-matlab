function a=read_array_nd(filename)
% a=rad_arrray_nd(filename)
% read nd matrix written by write_array_nd
d=load(filename);
if d(1) ~= 981,
    error('not written by write_arra_nd')
end
n=d(2);     % number of dimensions
bs=3;       % beginning of size
es=bs+n-1;  % index to to end of size
s=d(bs:es); % array size
if length(d)~=prod(s) + es
    error('wrong number of terms in the file')
end
a=reshape(d(es+1:end),s');
end


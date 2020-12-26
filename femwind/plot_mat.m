function plot_mat(X1,X2,M,j)
% plot_mat(X1,X2,M,i1,i2)
% in:
%   X1, X2  coordinates, cell size 3
%   M       matrix
%   j       index triple of the column of M to display
% example: plot_mat(X,XC,P,[5,5,1])

m=size(X1{1});
mm=prod(m);
n=size(X2{2});
nn=prod(n);
if any(size(M)~=[mm,nn])
    error('wrong matrix M size')
end
if ~isvector(j) | length(j) ~= 3,
    error('j must be vector length 3')
end
XX2=zeros(3,1);
for k=1:3, 
    XX2(k,1)=[X2{k}(j(1),j(2),j(3))];
end
jx = sub2ind(n,j(1),j(2),j(3));
ii=find(M(:,jx));
mm=length(ii);
XX1=zeros(3,mm);
MM=zeros(mm,1);
for i=1:mm
    ix = ii(i);
    [i1,i2,i3]=ind2sub(m,ix);
    for k=1:3, 
        XX1(k,i)=[X1{k}(i1,i2,i3)];
    end
    MM(i,1)=M(ix,jx);
end
hexa_plot(XX1,XX2,MM)

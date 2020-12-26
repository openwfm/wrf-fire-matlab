function hexa_plot(X1,X2,K)
% hexa_plot(X,K)
% plot matrix K on nodes X
% in
%   X1  size (3,m), node coordinates
%   X2  size (3,n), node coordinates
%   K  size (m,n), matrix
h = ishold;
[~,m]=size(X1);
[~,n]=size(X2);
if any([m,n]~=size(K))
    error('wrong sizes')
end
ml=8;
tol=1e-6;
symm = false;
if m==n,
    symm = norm(K-K','fro')<tol;
end
for j=1:n
    % plot3(X(1,j),X(2,j),X(3,j),'o','Linewidth',ml); hold on
    for i=1:m
        for k=1:3,xij{k}=[X1(k,i),X2(k,j)];end
        if symm,
            s = K(i,j)/sqrt(K(i,i)*K(j,j));
        else
            s = K(i,j)/norm([K(:,j);K(i,:)'],inf);
        end
        if s>tol,
            plot3(xij{:},'--k','Linewidth',s*ml); hold on
        elseif s<-tol
            plot3(xij{:},'--r','Linewidth',-s*ml); hold on
        end
    end
end
if ~h, hold off, end
grid on
end

function hexa_plot(X,K)
% hexa_plot(X,K)
% in
%   X  size (3,n), node coordinates
%   K  size (n,n), symmetric matrix
h = ishold;
[m,n]=size(X);
ml=10/max(abs(K(:)));
tol=1e-6;
for i=1:n
    plot3(X(1,i),X(2,i),X(3,i),'o','Linewidth',K(i,i)*ml);hold on
    for j=i+1:n
        for k=1:3,xij{k}=[X(k,i),X(k,j)];end
        if K(i,j)>tol,
            plot3(xij{:},'--k','Linewidth',K(i,j)*ml);
        elseif K(i,j)<-tol
            plot3(xij{:},'--r','Linewidth',-K(i,j)*ml);
        end
    end
end
if ~h, hold off, end
grid on
end

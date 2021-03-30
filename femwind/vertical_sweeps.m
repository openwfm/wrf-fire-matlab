function [x]=vertical_sweeps(K,F,X,x)

n = size(X{1});

for rb1=1:2
    for rb2=1:2
        for i1=rb1:2:n(1)
            for i2=rb2:2:n(2)
            % solving horizontal location i1 i2 and vertical line
                for i3=1:n(3)
                    ix = sub2ind(n,i1,i2,i3); 
                    x(ix) = x(ix) - K(ix,ix)\(K(:,ix)'*x - F(ix));
                    x(ix)
                end
            end
        end
    end
end
end

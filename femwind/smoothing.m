function x=smoothing(K,F,X,x,params)
% x=smoothing(K,F,X,x)
% one smoothing iteration
    n = size(X{1}); 
    it_type=sprintf('level %g smoothing',params.levels);
    switch params.smoothing
        case {'horizontal planes'}
            % disp('red-black vertical, planes horizontal')
            for rb3=1:2
                for i3=rb3:2:n(3)
                        [planex,planey]=ndgrid(1:n(1),1:n(2));
                        planez = i3*ones(n(1),n(2));
                        % solving horizontal layer
                        ix = sub2ind(n,planex,planey,planez); 
                        x(ix) = x(ix) - K(ix,ix)\(K(:,ix)'*x - F(ix));
                end
            end
        case {'vertical lines'}
            % disp('red-black relaxation horizontal, lines vertical')
            for rb1=1:2
                for rb2=1:2
                    for i1=rb1:2:n(1)
                        for i2=rb2:2:n(2)
                            % solving horizontal location i1 i2 and vertical line
                            ix = sub2ind(n,i1*onez,i2*onez,colz); 
                            x(ix) = x(ix) - K(ix,ix)\(K(:,ix)'*x - F(ix));
                        end
                    end
                end
            end
        case {'vertical sweeps'}
            % disp('red-black relaxation horizontal, down to up sweep vertical')
            x = vertical_sweeps(K,F,X,x);
        case '3D red-black'
            for rb1=1:2
                for rb2=1:2
                    for rb3=1:2    
                        for i1=rb1:2:n(1)
                            for i2=rb2:2:n(2)
                                for i3=rb3:2:n(3)
                                    ix = sub2ind(n,i1,i2,i3); 
                                    x(ix) = x(ix) - K(ix,ix)\(K(:,ix)'*x - F(ix));
                                end
                            end
                        end
                    end
                end
            end
        otherwise
            error(['smoothing ',params.smoothing,' unknown'])
    end
end
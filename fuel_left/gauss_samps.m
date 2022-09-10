function s = gauss_samps(n,split,split2)

s = zeros(5,n);
for i = 1:n
    if rand < split
        %faster ROS
        s(1:4,i) = 1000*randn(4,1);
    else
        %slower ros
        s(1:4,i) = 5000*randn(4,1);
    end
end
%set fuel time constants
for i = 1:n
    if rand < split2
        %slower burning fuels
        s(5,i) =  8 + 2000*rand;
    else
        %fast burning fuels
        s(5,i) =  8 + 10*rand;
    end
end


% figure,histogram(s(1,:))
% title('t00')


end % function

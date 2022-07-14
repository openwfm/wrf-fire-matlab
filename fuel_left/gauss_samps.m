
function s = gauss_samps(n,split)

s = zeros(5,n);
for i = 1:n
    if rand < split
        s(1:4,i) = 1000*randn(4,1);
    else
        s(1:4,i) = 5000*randn(4,1);
    end
end
%set fuel time constants
s(5,:) = 8 + 2000*rand(1,n);

% figure,histogram(s(1,:))
% title('t00')


end % function
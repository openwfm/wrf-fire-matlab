B=zeros(9);
for i=1:3:7
    B(i:i+2,i:i+2)=rand(3);
end
D = [1 1 1 0 0 0 0 0 0 
     0 0 0 1 1 1 0 0 0
     0 0 0 0 0 0 1 1 1];
% D = rand(2,4
C = [1 0 0 1 0 0 0 0 0
     0 0 0 1 0 0 0 1 0    
     0 0 0 0 0 0 1 0 0]
v0 = rand(9,1)
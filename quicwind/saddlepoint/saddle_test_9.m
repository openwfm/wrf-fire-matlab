B=zeros(6);
for i=1:2:6
    B(i:i+1,i:i+1)=rand(2);
end
D = [1 1 0 0 0 0 
     0 0 1 1 0 0
     0 0 0 0 1 1];
C = [1 0 1 0 0 0
     0 0 0 1 1 0    
     0 0 0 0 0 1]
v0 = rand(6,1)

saddle_sparse
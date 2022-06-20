%script creates data and trains the network. It then creates new data and
%tests the network

n = input_num('How many samples to create?',10000);
%create samples
s = time_samps(n);
%s(:,1:4) = rotate_cell(s(:,1:4));

%perform integrals
fprintf('Computing Integrals on data, be patient...\n')
%r = fuel_int(s);
r = fuel_quad(s,50);

fprintf('Now training the network \n')
train_simple

%create new set of data
s2 = time_samps(n);
tic 
r2 = fuel_quad(s2,50);
toc

%use network to evaulate integrals
tic
y2 = net(s2');
%set outlier values to zero or one
y2(y2>1)=1;
y2(y2<0)=0;
y2 = y2';
toc

%scatter plot of quadrature and network outputs
figure,scatter(r2,y2)
title('Network and integral evaluations')
xlabel('Integral values')
ylabel('Network values')

%histogram of differences
d = r2-y2;
figure,histogram(d)
title('Integral - Network')
xlabel('int - net')
ylabel('Instances')

fprintf('Maximum absolute difference is %f \n',max(abs(d)))
fprintf('Mean of difference is %f \n',mean(d))
fprintf('Variance of difference is %f \n',var(d))

    




function net_int = train_2_layers(s,r)
%function takes sample data an trains a neural network to aaproximate the
%integral used to compute the fuel fraction.
%inputs -
%   s -- samples from gauss_samps(1000,0.8), etc
%   r -- responses from fuel_quad(s,50), etc

%create CNN with two hidden layers of size 10
net_int = feedforwardnet([10 10])

%sigmoid activation function is used by default, but may be changed
%poslin --> ReLu activation function
%net_int.layers{1}.transferFcn = 'poslin';
%net_int.layers{2}.transferFcn = 'poslin';

%train the network on the data
[net_int, train_int] = train(net_int,s,r);


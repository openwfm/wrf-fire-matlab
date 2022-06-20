function s = fuel_samps(n)
%creates simulated data for testing and training the neural network to be
%used in module_fr_sfire_core.F
%Inputs -
%   n - number of samples to create.
%Output
%   s - nx5 matrix with rows representing times at the corners of a cell and
%   the fuel constant _fuel_cell_time. Each row is arranged as
%   [t00, t01 t10 t11 fuel_time_cell], where t00,...,t11 are taken to be
%   the time since ignition arranged as:
%       t01 ----- t11
%        |         |
%        |         |
%       t00 ----- t10

%matrix for [t00, t01 t10 t11 fuel_time_cell]
s = zeros(n,5);

%create data
for i = 1:n
    

function fq = fuel_int(s)
%computes fuel left using Matlab integral2 function
%
%
%
%                                  /\
%                                  |              
%                         fq =     | 1 -  exp(-T(x,y)/fuel_time_cell)) dxdy
%                                  |            
%                                 \/
%                                0<x<1
%                                0<y<1
%                             L(x,y)<=0
%Inputs:
%   s - nx5 matrix with rows 
%        [tign00 tign01 tign10 tign11 fuel_time_cell], where t00, t01, t10,
%        t11 are the time since ignition.
%Outputs:
%   fq - nx1 vector with fuel left

[n,~] = size(s);
fq = zeros(n,1);
%loop
for i = 1:n
    
    tign00 = s(i,1);
    tign01 = s(i,2);
    tign10 = s(i,3);
    tign11 = s(i,4);
    ftc = s(i,5);
    %don't integrate when cell has no fire in it.
    if any([tign00 tign01 tign10 tign11] > 0)
      %integrand made from interpolation of fire time,
      intt = @(x,y) 1 - ... 
                    exp( ...
                        -(((1.0-x).*( (1.0-y).*tign00 + y*tign01 ) + x.*( (1.0-y).*tign10 + y.*tign11 )) ...
                    + abs(((1.0-x).*( (1.0-y).*tign00 + y*tign01 ) + x.*( (1.0-y).*tign10 + y.*tign11 )))) ...
                     /2/ftc);
      
      %...
               %  +abs((1.0-x).*( (1.0-y).*tign00 + y*tign01 ) + x.*( (1.0-y).*tign10 + y.*tign11 )))/2/ftc;
    
      fq(i) = integral2(intt,0,1,0,1);
    end % if
end % for loop
end % function
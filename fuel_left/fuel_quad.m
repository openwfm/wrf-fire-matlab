function fq = fuel_quad(samp,m,n)
%computes quadrature
%samp - matrix with [tign00 tign01 tign10 tign11 time_now fuel_time_cell]
%or sp -- , mtrix with [tign00 tign01 tign10 tign11 fuel_time_cell]where
%   t00 is processed so that t00 = time_now - t00
%m - size of quadrature grid (m x m)
%n - number of integrals to evaluate

proc = 1; %
%proc = input_num('Is data processed? 1 = yes',0);

%transpose if samplles are from column vectors
sz = size(samp);
if sz(1)<sz(2)
    samp = samp';
end


%corners of cell [a,b]X[a,b]
a = 0.0;
b = 1.0;
%grid spacing
dx = (b-a)/m;
dy = dx;

if ~exist('n','var')
    n = length(samp);
end
fq = zeros(1,n);
%loops
tic
for o = 1:n % outer loop
    %fire parameters
    tign00 = samp(o,1);
    tign01 = samp(o,2);
    tign10 = samp(o,3);
    tign11 = samp(o,4);
    if proc == 1
        fuel_time_cell = samp(o,5);
    else
        time_now  = samp(o,5);
        fuel_time_cell = samp(o,6);
    end
    %integral sum --> s
    s = 0;
    %midpoint method
    for i = 1:m
        for j = 1:m
            x = a + dx/2.0 + (i-1.0)*dx;
            y = a + dy/2.0 + (j-1.0)*dy;
            % interpolate tign
            intt = (1.0-x)*( (1.0-y)*tign00 + y*tign01 ) + x*( (1.0-y)*tign10 + y*tign11 );
            % t is the time since fire arrived = time_now - intt
            if proc == 1
                t = intt;
            else
                t = time_now - intt;
            end
            if t >= 0
                bf = (1 - exp(-t/fuel_time_cell));
                % make sum update
                s = s + bf;
            end %if
        end % for j
    end % for i
    %multiply differencials
    fq(1,o) = s*dx*dy;

    %matlab integral2 method
    %intt = (1.0-x)*( (1.0-y)*tign00 + y*tign01 ) + x*( (1.0-y)*tign10 + y*tign11 );


end % for o
t = toc;
%fprintf('Elepsed time is %f , Each quad was %f \n', toc, toc/n);
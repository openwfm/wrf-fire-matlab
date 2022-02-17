function rv = ros_var(d,c,l,sig)
%gives standard deviation of estimate ros between two detections
%d - distance in meters
%c - time between detections in hours
%l - optional time of uniform disributiuon
%sig - standard deviation of geolocation error
if ~exist('l','var')
   l = 6;
   sig = 375;
end
if c < l
    rv = NaN;
else
    rv =  sqrt((d^2+2*sig^2)/l^2*log(c^2/(c^2-l^2)) - d^2/l^4*(log(((c+l)/(c-l))^l*((c^2-l^2)/(c^2))^c))^2)/3600;
end

end %function

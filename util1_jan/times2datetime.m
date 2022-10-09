function d=times2datetime(times)

if ~isvector(times),
    error('input must be vector');
end

timec = char(times(:)');

d = datetime(times,'InputFormat','yyyy-MM-dd_HH:mm:ss');

end
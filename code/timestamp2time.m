function t = timestamp2time(str)
% dt = textscan('2024-03-27 05:05:28', '%4d-%2d-%2d %2d:%2d:%2d');
dt = textscan(str, '%4d-%2d-%2d %2d:%2d:%2d');
t = datenum(datetime(dt{:}));
end
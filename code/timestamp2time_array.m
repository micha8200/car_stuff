function t = timestamp2time_array(chArray)
%%          
%           1234567890
% chArray = ['2024-03-27 05:05:28';
%            '2024-03-27 07:05:32'];
% uint8('0123456789') =
% 48   49   50   51   52   53   54   55   56   57
% s = chArray(:, 19)-48 + 10*(chArray(:, 18)-48);
% m = chArray(:, 16)-48 + 10*(chArray(:, 15)-48);
% dt = textscan('2024-03-27 05:05:28', '%4d-%2d-%2d %2d:%2d:%2d');

if iscellstr(chArray)
    chArray = char(chArray);
end
% will change a lot
s = chArray(:, 19) + 10*chArray(:, 18) -528;
m = chArray(:, 16) + 10*chArray(:, 15) -528;
h = chArray(:, 13) + 10*chArray(:, 12) -528;
D = chArray(:, 10) + 10*chArray(:, 9) -528;

% will not change at all during convertions, so take only first, the rest
% by diff!!!
M = chArray(:, 7) + 10*chArray(:, 6) -528;
Y = chArray(:, 4) + 10*chArray(:, 3) + 100*chArray(:, 2) + 1000*chArray(:, 1)-53328;

% numer of days since calendar t0 January 0, 0000
t0 = datenum(datetime(Y(1), M(1), D(1), h(1), m(1), s(1)));

t = t0 + D - D(1) + (h - h(1))/24 + (m - m(1))/24/60 + (s - s(1))/24/60/60;

end

function i = c2i(c)
i=c-48;
end
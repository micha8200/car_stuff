
RH = 0.8;
T = 35;


Rd = 287.05; % j/(kg*K)
Rv = 461.495; % j/(kg*K)
pd = RH*6.1094-exp(17.625*T/(T+243.04));% partial pressure dry air
pv = 101325-1909.4;% partial pressure water vapor


rho = (pd/Rd+pv/Rv)/(273+T)

% Temperature (°C)	Air Density (kg/m³)
Ts       = [10:5:40];
rhos    = [1.24210000000000;1.21880000000000;1.19570000000000;1.17280000000000;1.14970000000000;1.12630000000000;1.10240000000000];



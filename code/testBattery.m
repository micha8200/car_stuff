hold on
% plot(tb.Timestamp_IST_, tb.Power_kW_)
y = smooth(tb.Power_kW_, 0.1, 'loess');
plot(tb.Timestamp_IST_, y)
grid on
yyaxis right
y = smooth(tb.Elevation_m_, 0.1, 'loess');
plot(tb.Timestamp_IST_, y)


%%

% tb = readtable('LRW3E7FS8PC835950-2024-01-21 18-22-14-2024-01-21 19-26-03.csv');

loesmooth = @(y) smooth(y, 0.1, 'loess');
figure()

batkWh1 = tb.BatteryLevel___/100*58;
batkWh2 = tb.BatteryRange_km_*137/1e3;
hold on

colors = lines(7);


% plot(tb.Timestamp_IST_, batkWh1, 'o', 'Color', colors(1,:), 'displayn', 'byPercent');
plot(tb.Timestamp_IST_, loesmooth(batkWh1), 'o', 'Color', colors(1,:), 'displayn', 'byPercent');
plot(tb.Timestamp_IST_, batkWh2, '.', 'Color', colors(2,:), 'displayn', 'byRange');
legend show
grid on
grid minor
ylim([0 60])
%%
figure()
kwhEnergyLevel = loesmooth(batkWh2);

t = tb.Timestamp_IST_;


x = t;
x = tb.Odometer_km_ - tb.Odometer_km_(1);

plot(x, kwhEnergyLevel, '-', 'displayn', 'BatteryState');

% t = datenum(tb.Timestamp_IST_, 'yyyy-mm-dd HH:MM:SS');


tsec    = datenum(t)*86400;
dt      = diff(tsec);
pow1    = tb.Power_kW_(2:end);
dE      = dt.*pow1*1000/60/60; % Wsec->Wh
Ecum    = cumsum(dE);

plot(x(2:end), Ecum/1e3, '-', 'displayn', 'EnergyLoss');
% ylim([0 60])
ylabel('kWh');

%% altitude, velocity 

figure
deltaDist       = diff(smooth(tb.Odometer_km_, 0.1, 'sgolay'));
EnPerDist       = smooth(dE./deltaDist, 0.1, 'sgolay');
deltaAlt        = diff(tb.Elevation_m_);
csum_delta_alt  = cumsum(abs(deltaAlt));
efficiency      = 0.9;
cume_loss_alt   = (1-efficiency) * mass*9.81*csum_delta_alt * kwh_per_joule; %wh
eloss_alt       = abs(deltaAlt).*(1-efficiency) * mass*9.81 * wh_per_joule;
EnPerDist_wo_altloss     = smooth((dE-sign(dE).*eloss_alt)./deltaDist, 0.1, 'sgolay');

mHeight = smooth(tb.Elevation_m_, 0.1, 'loess');
plot(x, mHeight, '-', 'displayn', 'altitude[m]');

yyaxis right
hold on
plot(x(2:end), EnPerDist, '-', 'displayn', 'Usage[Wh/km]');
plot(x(2:end), EnPerDist_wo_altloss, '-', 'displayn', 'delta alt loss');
% ylim([0 60])
ylabel('Wh/km');


% [deltaAlt round(dE, 1) round(sign(dE).*eloss_alt, 1)]
legend show


%% calculate mechanical  + battery energy

file = 2;

save_close = 0;

if file==1 % 86km
    tb = readtable('LRW3E7FS8PC835950-2024-01-21 18-22-14-2024-01-21 19-26-03.csv'); 
elseif file==2 % 152km
    tb = readtable('C:\Users\MICHACHA\Downloads\LRW3E7FS8PC835950-2024-01-17 16-45-14-2024-01-17 18-47-07.csv');
elseif file==3 % 184km
    tb = readtable('C:\Users\MICHACHA\Downloads\LRW3E7FS8PC835950-2023-12-09 16-09-07-2023-12-09 17-58-28.csv'); % tamra ~185km
elseif file==4 % 125km
    tb = readtable('C:\Users\MICHACHA\Downloads\LRW3E7FS8PC835950-2023-12-09 08-39-00-2023-12-09 09-49-31.csv'); % haifa 
end


mass            = 1750+90;
kwh_per_joule   = 2.77778e-7;
wh_per_joule    = 2.77778e-4;
nstd            = 2; % 97% of population
epa_wh_per_km   = 137; % 60/437

% vel_bins = [-inf 10:10:120 inf];
% vel_bins = [0 20:20:120 inf];
% vel_bins = [0 20 50 70 90 100 110 120 125 130 inf];
% vel_bins = [0 40:20:100 110:10:130 inf];
vel_bins = [0 40:10:130 inf];

vel_bins_mid    = [(vel_bins(1:end-2) + vel_bins(2:end-1))/2 nan];
mDist           = tb.Odometer_km_ - tb.Odometer_km_(1);
batkWh2         = tb.BatteryRange_km_*epa_wh_per_km/1e3;
kwhEnergyLevel  = loesmooth(batkWh2);
mpsVel          = tb.Speed_km_h_/3.6; %smooth(tb.Speed_km_h_/3.6, 0.05, 'loess');
mAlt            = smooth(tb.Elevation_m_, 0.05, 'sgolay');
E_k             = 0.5*mass*mpsVel.^2; %j
E_m             = mass*9.81*mAlt; % j
E_bat           = kwhEnergyLevel; % kWh
E_total         = smooth(E_k*kwh_per_joule + E_m*kwh_per_joule + E_bat, 0.1, 'loess');
csum_delta_alt  = cumsum(abs(diff(tb.Elevation_m_)));
efficiency      = 0.9;
cume_loss_alt   = (1-efficiency) * mass*9.81*csum_delta_alt * kwh_per_joule; %wh

% -------------- plot energy use --------------
figure()
hold on
% plot(mDist(2:end), E_total(2:end)+cume_loss_alt, 'displayname', 'E_{kin}+E_{pot}+E_{battery}[kWh]+E_{regenloss}')
plot(mDist, E_total, 'displayname', 'E_{kin}+E_{pot}+E_{battery}[kWh]')
plot(mDist, E_bat, 'displayname', 'E_{battery}[kWh]')
% plot(mDist, E_k*kwh_per_joule + E_m*kwh_per_joule, 'displayname', 'E_{mechanizal}[kWh]')

ylabel('Energy [kWh]')
legend show

yyaxis right

plot(mDist, mpsVel*3.6, '.-', 'displayname', 'velocity[km/h]')
legend show
ylabel('Velocity [km/h]')


ix = discretize(tb.Speed_km_h_, vel_bins);

cix = accumarray(ix, [1:length(ix)]', [length(vel_bins)-1 1], @(x) {x});

dE_dx = gradient(E_total)./gradient(mDist);

vels = zeros(1, length(cix))';
vels_std = zeros(1, length(cix))';

slope = zeros(1, length(cix))';
slope_std = zeros(1, length(cix))';

bin_num = [1:length(cix)]';

for i=1:length(cix)
    ix = cix{i};
    
    vels(i) = mean(mpsVel(ix)*3.6);
    vels_std(i) = std(mpsVel(ix)*3.6);
    
    temp = dE_dx(ix);
    temp = temp(~isinf(temp));
    slope(i) = -mean(temp)*1e3;
    slope_std(i) = std(temp)*1e3;
    
end

tb2 = table();
tb2.bin = bin_num;
tb2.n = cellfun(@length, cix);
tb2.Vels = char([num2str(vel_bins(1:end-1)') repmat(' - ', [length(cix) 1]) num2str(vel_bins(2:end)')]);
vels = round(vels);
ncoeff = nstd * sqrt(tb2.n);

slope = round(slope);
vels_std = ceil(vels_std);
slope_e = ceil(2*slope_std./ncoeff);

sl1 = slope-slope_e;
sl2 = slope+slope_e;
tb2.EPD = char([num2str(sl1) repmat(' - ', [length(cix) 1]) num2str(sl2)]);
tb2.nominalEPD = num2str(slope);
sprintf('%1.0f-%1.0f\n', vel_bins(1:end-1), vel_bins(2:end))

disp(tb2)

% -------------- plot efficiency --------------
f = figure;
hold on
ok = (slope_e./abs(slope))<0.4;
errorbar(vel_bins_mid(ok), slope(ok), slope_e(ok))
plot([10:10:200], calcEnLoss(mass, [10:10:200]), 'displayn', 'Theoretical')


ylabel('Energy Usage [Wh/km]')
xlabel('Speed [km/h]')

ylim([100 200]);
grid minor

yyaxis right

b = bar(vel_bins_mid, tb2.n/sum(tb2.n)*100, 'displayn', 'speed', 'EdgeAlpha', 0.1, 'FaceAlpha', 0.1);
ylim([0 80]);
ylabel('Occurance [%]')
grid on
box on

filename = sprintf('%s_%1.0fkm.png', datestr(tb.Timestamp_IST_(1), 'yyyy_mm_dd_HHMM'), range(tb.Odometer_km_))

insideTempFld = find(contains(lower(tb.Properties.VariableNames), 'inside'));
outsideTempFld = find(contains(lower(tb.Properties.VariableNames), 'outside'));
inTemp = mean(tb.(tb.Properties.VariableNames{insideTempFld}));
outTemp = mean(tb.(tb.Properties.VariableNames{outsideTempFld}));

title(sprintf('Economy\nInside:  %1.1fC\nOutside: %1.1fC', inTemp, outTemp))

if save_close
    saveas(f, filename);

    close all
end

%% cumulative delta_alt


sum_delta_alt = sum(abs(diff(tb.Elevation_m_)));

csum_delta_alt = cumsum(abs(diff(tb.Elevation_m_)));
efficiency = 0.9;

e_loss_total = (1-efficiency) * mass*9.81*sum_delta_alt * kwh_per_joule; %wh
fprintf('cumulative delta alt = %1.0f => loss= %1.1fkWh\n', sum_delta_alt, e_loss_total)



%% theoretical usage limit

% calcEnLoss(mass, [100:10:130])

%%

tbbat = readtable('C:\Users\MICHACHA\Downloads\LRW3E7FS8PC835950-battery.csv', 'PreserveVariableNames', 0);

ix2rem = [10 22:24 29 30 31:35 41 43];

% tbbat(ix2rem,:) = [];

kWh_bat_rng = tbbat.MaxRange_km_*epa_wh_per_km/1e3;
kWh_bat_bat = tbbat.UsableCapacity_kWh_;
odo = tbbat.Odometer_km_;
ok              = true(size(odo));
ok(ix2rem)      = 0;

figure();
hold on
% plot(tbbat.Odometer_km_, kWh_bat_rng, ':', 'displayn', 'from\_range')
% plot(tbbat.Odometer_km_, smooth(kWh_bat_rng, 0.1), 'displayn', 'from\_range')
plot(tbbat.Odometer_km_, kWh_bat_bat, ':', 'displayn', 'from\_cap')
plot(tbbat.Odometer_km_(ok), kWh_bat_bat(ok), 'ok', 'displayn', 'from\_cap')
plot(tbbat.Odometer_km_(ok), smooth(kWh_bat_bat(ok), 0.2), 'displayn', 'from\_cap\_smooth')

% polyfit nth order

n = 1;
p = polyfit(tbbat.Odometer_km_(ok), kWh_bat_bat(ok), n);

%y_p

% x = tbbat.Odometer_km_(ok)';
x = 0:1e3:150e3;
% y = fliplr(p)*transpose(double(x).^[0:n]);
y = double(x').^[0:n]*p(end:-1:1)';

plot(x, y, 'displayn', 'from\_cap\_polyfit')

legend show

infusefig

%%

% N = 10;
clc
for N = [5 10 50 100 1000 10000 100000]
    vals = randn([1 N]);
    fprintf('N=%1.0f: %1.2f ± %1.2f\n', N, mean(1+vals), 2*std(vals)/sqrt(N)); % /sqrt(N)
end


%%
function en_loss = calcEnLoss(mass, vel)
% mass[kg]
% vel[km/h]
rho = 1.21; % kg/m^3
% rho = 1.204;
drag_coefficient = 0.23;
front_cross_section = 2.3;
rolling_resistance = 0.011; % tires+road

% drag
f_drag = front_cross_section*drag_coefficient*0.5*rho*(vel/3.6).^2; % N
enloss_drag = f_drag*100000/1000/3600*10; % Wh/km

% roll
enloss_roll = mass*9.81*rolling_resistance*100000/1000/3600*10; % Wh/km

en_loss = round(enloss_roll + enloss_drag);

end
function plotDrive(csvFile, varargin)

% get drive/state table from python script:
% C:\Users\micha\Documents\Python\Tesla\getStatesTocsv.py

% ideas: write each drive as a clean slate, so that it will become an
% energy use segment. Then enable cumulative addition of such segments, to
% allow multiple day-charge graphs, and statistics!!


%
if nargin<2
    issave = false;
    isclose = false;
    plotrawecon = false;
    isclearb4 = true;
else
    issave = any(contains(varargin, '-save'));
    isclose = any(contains(varargin, '-close'));
    isclearb4 = any(contains(varargin, '-clear'));
    plotrawecon = any(contains(varargin, '-raw'));
end

if nargin<1
    [csvFile, path]=uigetfile('C:\Users\micha\Documents\Python\Tesla\*.csv');
%     [csvFile, path]   = uigetfile('*.csv,*.txt');
    csvFile = fullfile(path, csvFile);
    if isempty(csvFile)
        return
    end
    fprintf('plotDrive %s\n', csvFile)
end
tb = readtable(csvFile); 
tb.Properties.VariableNames = strrep(tb.Properties.VariableNames, '_', '');

if isclearb4
    close(findobj(0, 'Name', 'driveEfficiency'));
    close(findobj(0, 'Name', 'driveEnergy'));
end

%%

% find idle energy loss:

ispark = ismember(lower(tb.ShiftState), 'p');

tbp = tb(ispark, :);
uid = unique(tbp.Odometerkm);
[~, ixu] = ismember(tbp.Odometerkm, uid);

deltaBatRanges = accumarray(ixu, tbp.BatteryRangekm, [], @(x) max(x) - min(x));

% ranges dropped by deltaBatRanges in odo states uid!!!
tb = tb(ismember(lower(tb.ShiftState), 'd'), :);
% optional: correct ALL odo states above each uid by + its' deltaBatRanges

for i=1:length(uid)
    ok = tb.Odometerkm>uid(i);
    tb.BatteryRangekm(ok) = tb.BatteryRangekm(ok) + deltaBatRanges(i);
end
% pratition into odometer bins


%%

tb = tb(ismember(lower(tb.ShiftState), 'd'), :);

assignin('base', 'tb', tb);
%% parameters

mass            = 1750+250;
% mass            = 1750;
loesmooth       = @(y) smooth(y, 0.1, 'loess');
EconomDispLims = [0 300]; % wh/km
econ_thresh     = 400; % wh/km threshold to calc (mainly to remove diverging data due to delta_Dist->0)
kmDist         	= tb.Odometerkm - tb.Odometerkm(1);
%% Energy

kwhFullBatteryTemp = max(tb.EnergyRemainingkWh);
batkWh1 = tb.BatteryLevel/100*kwhFullBatteryTemp;
E_real_r = tb.BatteryRangekm*GP.epa_wh_per_km/1e3;
kwhEnergyLevel = loesmooth(E_real_r);

if ismember('TimestampUTC', tb.Properties.VariableNames)
    t       = tb.TimestampUTC; % days? 
    tsec    = datenum(t)*86400;
elseif ismember('TimestampIDT', tb.Properties.VariableNames)
%     tsec = datenum('2025-08-11 06:14:51', 'yyyy-mm-dd HH:MM:SS');
    t0 = second(tb.TimestampIDT);
    t1 = minute(tb.TimestampIDT);
    t2 = hour(tb.TimestampIDT);
    t3 = day(tb.TimestampIDT);
    
    tsec = t0 + 60*t1 + 60*60*t2 + 60*60*24*t3;
    
%     tsec = zeros(height(tb));
%     for i=1:length(tsec)
%         tsec(i) = datenum(tb.TimestampIDT(i), 'yyyy-mm-dd HH:MM:SS');
%     end
end


dt      = diff(tsec);
pow1    = tb.PowerkW(2:end);
dE      = dt.*pow1*1000/60/60; % Wsec->Wh

Ecum    = cumsum(dE);



%% altitude, velocity 

f(1) = figure;
f(1).Name = 'driveEfficiency';

deltaDist       = diff(smooth(tb.Odometerkm, 0.1, 'sgolay'));
wh_km_econ      = dE./deltaDist;
wh_km_econ(wh_km_econ<-econ_thresh|wh_km_econ>econ_thresh) = nan;
EnPerDist       = smooth(wh_km_econ, 0.1, 'sgolay');
EnPerDistS      = smooth(wh_km_econ, 0.5, 'lowess');
deltaAlt        = diff(tb.Elevationm);
csum_delta_alt  = cumsum(abs(deltaAlt));
efficiency      = 0.9;
cume_loss_alt   = (1-efficiency) * mass*9.81*csum_delta_alt * GP.kwh_per_joule; %wh
eloss_alt       = abs(deltaAlt).*(1-efficiency) * mass*9.81 * GP.wh_per_joule;
EnPerDist_wo_altloss     = smooth((dE-sign(dE).*eloss_alt)./deltaDist, 0.1, 'sgolay');

mHeight         = tb.Elevationm;
mHeight_sm      = smooth(mHeight, 0.1, 'loess');
plot(kmDist, mHeight, '-', 'displayn', 'altitude[m]');
xlabel('DriveRange[km]')
ylabel('Altitude[m]')

yyaxis right
hold on
if plotrawecon
    plot(kmDist(2:end), EnPerDist, '-', 'displayn', 'Econ[Wh/km]');
end
plot(kmDist(2:end), EnPerDistS, ':', 'displayn', 'Econ(smooth)[Wh/km]');
colormap lines(3)
ylim(EconomDispLims)
% plot(x(2:end), EnPerDist_wo_altloss, '-', 'displayn', 'delta alt loss');
% ylim([0 60])
ylabel('Wh/km');
grid on
grid minor

% [deltaAlt round(dE, 1) round(sign(dE).*eloss_alt, 1)]
legend show

%%

f(2)= figure();
f(2).Name = 'driveEnergy';


hold on
box on
grid on
grid minor

if ~any(ismember(tb.Properties.VariableNames, 'Speedkmperh')) && any(contains(lower(tb.Properties.VariableNames), 'speed'))
    ix = find(contains(lower(tb.Properties.VariableNames), 'speed') & contains(lower(tb.Properties.VariableNames), 'km'));
    tb.Properties.VariableNames{ix} = 'Speedkmperh'; %#ok<FNDSB>
end
    
E_real_r         = tb.BatteryRangekm*GP.epa_wh_per_km/1e3;
kwhEnergyLevel  = loesmooth(E_real_r);
kwhE0           = E_real_r(1);
E_epa           = kwhE0 - GP.epa_wh_per_km*kmDist/1e3; % including height
kwhDeltaEnByHeight = mass*9.81*(mHeight - mHeight(1)) * GP.kwh_per_joule;
kwhCurrentSpeedEnergy = 0.5*mass*(tb.Speedkmperh/3.6).^2*GP.kwh_per_joule;

% current height corrected EPA drive Energy state
E_EPA_HC        = E_epa-kwhDeltaEnByHeight;
E_REAL_HC       = E_real_r + kwhDeltaEnByHeight;

% current height and current velocity corrected EPA Energy drive state
E_EPA_HC_VC     = E_epa-kwhDeltaEnByHeight-kwhCurrentSpeedEnergy;
E_REAL_HC_VC     = E_real_r+kwhDeltaEnByHeight+kwhCurrentSpeedEnergy;

plot(kmDist, E_epa, ':k', 'LineWidth', 2, 'DisplayName', 'E_{EPA} flat')
plot(kmDist, E_EPA_HC, ':', 'LineWidth', 2, 'DisplayName', 'E_{EPA} w Epot')
% plot(kmDist, E_EPA_HC_VC, ':', 'LineWidth', 2, 'DisplayName', 'E_{EPA} w Epot+Ekin')
plot(kmDist, E_real_r, '-', 'LineWidth', 2, 'DisplayName', 'E_{battery}')
% plot(kmDist, E_REAL_HC, '-', 'LineWidth', 2, 'DisplayName', 'E_{battery} w Epot')
p = plot(kmDist, E_REAL_HC_VC, '-', 'LineWidth', 2, 'DisplayName', 'E_{battery} w Epot+Ekin');
p.UserData = tb;


title(sprintf('Energy state\ndifference from EPA:%1.2fkWh in %1.0fkm\nTotal efficiency: %1.0f Wh/km', E_EPA_HC_VC(end)-E_real_r(end), range(kmDist), range(E_real_r*1e3)/range(kmDist)))
ylabel('Energy [kWh]')
xlabel('Distace traveled [km]')
lg = legend('show');
lg.ItemHitFcn = @LegendItemHitFnc;

yyaxis right

plot(kmDist, tb.Speedkmperh, '.-', 'displayname', 'velocity[km/h]')
legend show
ylabel('Velocity [km/h]')

legend show

%%

if issave
    saveas(f, filename);
end
if isclose
    close all
end


end


%% service functions

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

function plotEnergy(csvFile, varargin)

% plotEnergy my_file.csv -pe -pd -idle3kw -temp35
% appendTables
if nargin<2
    issave = false;
    isclose = false;
    plotDrive = false;
    plotEnergy = false;
else
    issave = any(contains(varargin, '-s'));
    isclose = any(contains(varargin, '-c'));
    plotDrive = any(contains(lower(varargin), '-plotd'));
    plotEnergy = any(contains(lower(varargin), '-plote'));
end
if any(contains(varargin, '-idle'))
    str = varargin{find(contains(varargin, '-idle'), 1, 'last')};
    str = strrep(str, '-idle', '');
    if contains(str, 'k')
        str = strrep(str, 'k', '');
        c = 1e3;
    else
        c = 1;
    end
    str = strrep(str, 'w', '');
    str = strrep(str, 'W', '');
    const_pow = str2num(str)*c;
else
    const_pow = [];
end
if any(contains(varargin, '-temp'))
    str = varargin{find(contains(varargin, '-temp'), 1, 'last')};
    str = strrep(str, '-temp', '');
    outTemp = str2num(str);
else
    outTemp = [];
end
if ~plotEnergy && ~plotDrive
    fprintf('plot not chosen! choose:\n-pd   : cumulative drives energy usage\n-pe   : cumulative data energy efficiency\n\n')
%     return
end
if nargin<1
    csvFile   = uigetfile('LRW3E7FS8PC835950*.csv');
    if isempty(csvFile)
        return
    end
end
% modelType = csvFile(4);
% modelYear = strfind('MNPRSTVWXY', csvFile(10)) + 2020;
% % 'MNPRSTVWXY' 'M' - 2021
[~, f, ~]=fileparts(csvFile);
vin = f(1:17);
mc = tesla_models_specs(vin); % model consants

tb = readtable(csvFile); 
tb.Properties.VariableNames = strrep(tb.Properties.VariableNames, '_', '');

% 
% if ~any(ismember(tb.Properties.VariableNames, 'Speedkmperh')) && any(contains(lower(tb.Properties.VariableNames), 'speed'))
%     ix = find(contains(lower(tb.Properties.VariableNames), 'speed') & contains(lower(tb.Properties.VariableNames), 'km'));
%     tb.Properties.VariableNames{ix} = 'Speedkmperh'; %#ok<FNDSB>
% end

%% change to segments

% segTime
% segDist = [0;];

% tb = readtable('MY_drives\XP7YGCEK5SB651296_merged.csv');

% t = tb.Timestamp_IDT_;
% t = double(second(t) + 60*minute(t) + 60*60*hour(t) + 60*60*24*day(t));
% [0; diff(t)>30];
% find lines to remove

% some tables have their power/energy NaN
ixnan = isnan(tb.EnergyRemainingkWh) & ~isnan(tb.BatteryRangekm);
tb.EnergyRemainingkWh(ixnan) = tb.BatteryRangekm(ixnan)/mc.epa_range*mc.kwh_in_new_battery;

tb.kwh_deltaE   = [0; diff(tb.EnergyRemainingkWh)];
tb.km_deltaX    = [0; diff(tb.Odometerkm)];
tb.s_deltaT     = [0; diff(tb.TimestampIDT)];

% ix = find(seconds(tb.s_deltaT)>30 | tb.km_deltaX>2 | tb.kwh_deltaE>2);
ix = find(seconds(tb.s_deltaT)>50 | tb.km_deltaX>5 | tb.kwh_deltaE>2);

fprintf('total paths before filter : %1.0fkm\n', sum(tb.km_deltaX));
tb(ix,:) = [];
fprintf('total paths after filter  : %1.0fkm\n', sum(tb.km_deltaX));

if any(contains(lower(varargin), 'vel'))
    if any(contains(varargin, '+'))
        tb = sortrows(tb, 'Speedkmh', 'descend');
    else
        tb = sortrows(tb, 'Speedkmh');
    end

end

%%

loesmooth       = @(y) smooth(y, 0.1, 'loess');
mass            = mc.kg_vehicle_mass;
if isempty(const_pow)
    const_pow       = 1500; % constant energy drain (power): AC, powertrain, inefficiencies. all independant of veocity
end
kwh_per_joule   = 2.77778e-7;
wh_per_joule    = 2.77778e-4;
nstd            = 2; % 97% of population
epa_wh_per_km   = mc.epa_wh_per_km; % 60/437

% vel_bins = [-inf 10:10:120 inf];
% vel_bins = [0 20:20:120 inf];
% vel_bins = [0 20 50 70 90 100 110 120 125 130 inf];
% vel_bins = [0 40:20:100 110:10:130 inf];
vel_bins = [0 20:10:140 inf];
vel_bins = [0:10:140 inf];
vel_bins = [-inf 15:10:145 inf];

vel_bins_mid    = [(vel_bins(1:end-2) + vel_bins(2:end-1))/2 nan];
% kmDist           = tb.Odometerkm - tb.Odometerkm(1);
% batkWh2         = tb.BatteryRangekm*epa_wh_per_km/1e3;

kmDist      	= cumsum(tb.km_deltaX);
batkWh2         = cumsum(tb.kwh_deltaE);

kwhEnergyLevel  = loesmooth(batkWh2);
mpsVel          = tb.Speedkmh/3.6; %smooth(tb.Speed_km_h_/3.6, 0.05, 'loess');
mAlt            = smooth(tb.Elevationm, 0.05, 'sgolay');
E_k             = 0.5*mass*mpsVel.^2; %j
E_m             = mass*9.81*mAlt; % j
E_bat           = kwhEnergyLevel; % kWh
E_epa100p       = - epa_wh_per_km*kmDist/1e3;
E_total         = smooth(E_k*kwh_per_joule + E_m*kwh_per_joule + E_bat, 0.1, 'loess');
csum_delta_alt  = cumsum(abs(diff(tb.Elevationm)));
efficiency      = 0.9;
cume_loss_alt   = (1-efficiency) * mass*9.81*csum_delta_alt * kwh_per_joule; %wh

%% -------------- plot energy use --------------

% plot1 = 0;
if plotDrive
    
    cc = lines(7);
    figure('name', mc.name)
    hold on
    fprintf('Energy state given:\n1.K+P+E\n2.BatteryState\n3.100%%efficiency(EPA)4.add option to change Yax to altitude?\n\n');
    
    % plot(mDist(2:end), E_total(2:end)+cume_loss_alt, 'displayname', 'E_{kin}+E_{pot}+E_{battery}[kWh]+E_{regenloss}')
    
    de = E_total-E_total(1);
    plot(kmDist, de, 'linewidth', 2, 'displayname', 'K+P+E_{battery} [kWh]', 'Color', cc(1,:))
    
    ok = true(size(kmDist));
%     ok = kmDist<175;
%     ok = kmDist>185;
%     ok = kmDist>68 & kmDist<166;        
    p = polyfit(kmDist(ok), de(ok), 1);
    p(2) = 0;
    x = [0 (-mc.kwh_in_new_battery -p(2))/p(1)];
    y = p(1)*x + p(2);
    plot(x, y, ':', 'linewidth', 1.2, 'displayname', 'mean usage', 'Color', cc(1,:))

    txt = sprintf('expected range %1.1fkWh for %1.0fkm', abs(y(2)), abs(x(2)));
    t = text(x(2), y(2), txt);
    t.VerticalAlignment = 'bottom';
    t.HorizontalAlignment = 'right';

%     wgtmean_speed = sum(tb.Speedkmh.*kmDist)./(sum(kmDist));
%     title(sprintf('TotalEnergy (speed, height invariant) vs time\n actual econ: %1.0f wh/km\naverege speed %1.0fkm/h', abs(p(1)*1e3), wgtmean_speed))
    
    pc80 = prctile(tb.Speedkmh, 80);
    title(sprintf('TotalEnergy (speed, height invariant) vs time\n actual econ: %1.0f wh/km\ntop 20%% speed %1.0fkm/h', abs(p(1)*1e3), pc80))
    
    
    % plot(kmDist, E_bat, 'linewidth', 2, 'displayname', 'E_{battery}[kWh]')
    x = [0:10:mc.epa_range mc.epa_range+eps];
    y = -x*mc.epa_wh_per_km/1e3;
    plot(x, y, ':k', 'linewidth', 1.5, 'displayname', 'EPA\_Efficiency')   
    
    x = [0:10:mc.wltp_range mc.wltp_range+eps];
    y = -x*mc.wltp_wh_per_km/1e3;
    plot(x, y, ':k', 'linewidth', 1.5, 'displayname', 'WLTP\_Efficiency')  
    
    
    plot([0 mc.epa_range mc.wltp_range], -mc.kwh_in_new_battery * [1 1 1], 'o:r', 'DisplayName', 'one drive full battery', 'LineWidth', 2)
    text(mc.epa_range, -mc.kwh_in_new_battery, 'EPA', 'VerticalAlignment', 'top')
    text(mc.wltp_range, -mc.kwh_in_new_battery, 'WLTP', 'VerticalAlignment', 'top')
    grid on
    grid minor
    % plot(mDist, E_k*kwh_per_joule + E_m*kwh_per_joule, 'displayname', 'E_{mechanizal}[kWh]')

    
    
    
    ylabel('Energy [kWh]')
    
    
    legend show

    yyaxis right
    rax = gca;
    rax.YColor = cc(3,:);
    
    plot(kmDist, mpsVel*3.6, '.-', 'displayname', 'velocity[km/h]', 'Color', cc(3,:))
    legend show
    ylabel('Velocity [km/h]')
    
    
    xlim([0 mc.wltp_range])

end
ix = discretize(tb.Speedkmh, vel_bins);

cix = accumarray(ix, [1:length(ix)]', [length(vel_bins)-1 1], @(x) {x});

dE_dx = gradient(E_total)./gradient(kmDist);

vels = zeros(1, length(cix))';
vels_std = zeros(1, length(cix))';

slope = zeros(1, length(cix))';
slope_std = zeros(1, length(cix))';

bin_num = [1:length(cix)]';

kmDeltaX = tb.km_deltaX;
kwhDeltaE = tb.kwh_deltaE;
kwhDeltaE = [0;diff(E_total)];
sDeltaT = tb.s_deltaT;

for i=1:length(cix)
    ix = cix{i};

    x = cumsum(kmDeltaX(ix));
    E = cumsum(kwhDeltaE(ix));
    
    if length(x)>4
        p = polyfit(x, E, 1);
    else
        fprintf('#%1.0g points in %1.0f-%1.0f  km/h interval\n', length(x), vel_bins(i), vel_bins(i+1))
        p = nan(1, 2);
    end
    slope(i) = -p(1)*1e3;
    
    pE = p(2) + p(1)*x;
    dE = pE - E;
    slope_std(i) = sqrt(sum(dE.^2/length(E)))*1e3;
%     vel_bin_std = 
%     slope_std(i) = std(temp)*1e3;
    
    
    vels(i) = mean(mpsVel(ix)*3.6);
    vels_std(i) = std(mpsVel(ix)*3.6);

    
%     temp = dE_dx(ix);
%     temp = temp(~isinf(temp));
%     slope(i) = -mean(temp)*1e3;
%     slope_std(i) = std(temp)*1e3;

end

tb2 = table();
tb2.bin = bin_num;
tb2.n = cellfun(@length, cix);
tb2.Vels = char([num2str(vel_bins(1:end-1)') repmat(' - ', [length(cix) 1]) num2str(vel_bins(2:end)')]);
vels = round(vels);
tb2.meanVel = vels;
ncoeff = nstd * sqrt(tb2.n);

slope = round(slope);
vels_std = ceil(vels_std);
slope_e = ceil(2*slope_std./ncoeff);

sl1 = slope-slope_e;
sl2 = slope+slope_e;
tb2.EPD = string([num2str(sl1) repmat(' - ', [length(cix) 1]) num2str(sl2)]);
tb2.nominalEPD = num2str(slope);

% sprintf('%1.0f-%1.0f\n', vel_bins(1:end-1), vel_bins(2:end))

disp(tb2)



%% -------------- plot efficiency --------------

if plotEnergy
f = figure('name', mc.name);
hold on
ok = ~isnan(slope_e) & slope>0;
% ok = (slope_e./abs(slope))<0.4 & vel_bins_mid'>20;

errorbar(vel_bins_mid(ok), slope(ok), 2*slope_e(ok), 'DisplayName', 'calculated efficiency-velocity');


vels2test = [10:10:200];

const_pow_temp = const_pow + diff(0.5*mass*(vels+vels_std.*[-1 1]).^2, [], 2)*wh_per_joule; 
const_pow2 = zeros(size(vels2test));
ix = discretize(vels, vels2test);
ok = ~isnan(ix);
const_pow2(ix(ok)) = const_pow_temp(ix(ok));

% energy loss for vel changes inside the bin
% higher velocity STD = highere regen/braking use in this speed category
if isempty(outTemp)
    outTemp = mean(tb.OutsideTempC);
end
[Eloss, Eloss_c] = calcEnLoss(mc, vels2test, const_pow, outTemp, mean(tb.InsideTempC));
plot(vels2test, Eloss, 'displayn', 'theoretical efficiency')
plot(vels2test, Eloss_c, 'displayn', sprintf('loss to const power (%1.1fkW)', const_pow/1e3))
[~, ixmin]=min(Eloss);
vel_epa     = interp1(Eloss(ixmin:end), vels2test(ixmin:end), mc.epa_wh_per_km);
vel_wltp    = interp1(Eloss(ixmin:end), vels2test(ixmin:end), mc.wltp_wh_per_km);
x = [vels2test([ 1 end]) nan vel_epa vel_epa];
y = [mc.epa_wh_per_km mc.epa_wh_per_km nan mc.epa_wh_per_km 0];
plot(x, y, ':', 'displayn', 'EPA', 'linewidth', 1.5)
x = [vels2test([ 1 end]) nan vel_wltp vel_wltp];
y = [mc.wltp_wh_per_km mc.wltp_wh_per_km nan mc.wltp_wh_per_km 0];
plot(x, y, ':', 'displayn', 'WLTP', 'linewidth', 1.5)
text(vel_epa, 0, sprintf('EPA-%1.0fkm/h', vel_epa), 'Rotation', 30);
text(vel_wltp, 0, sprintf('WLTP-%1.0fkm/h', vel_wltp), 'Rotation', 30);
ss1 = scatter(tb.Speedkmh, -dE_dx*1e3, '.', 'DisplayName', 'efficiency-velocity samples');

ylabel('Energy Usage [Wh/km]')
xlabel('Speed [km/h]')

ylim([0 400]);
xlim([0 200]);
grid minor

yyaxis right

b = bar(vel_bins_mid, tb2.n/sum(tb2.n)*100, 'displayn', 'speed bins', 'EdgeAlpha', 0.1, 'FaceAlpha', 0.1);
ylim([0 80]);
ylabel('Occurance [%]')
grid on
box on





insideTempFld   = find(contains(lower(tb.Properties.VariableNames), 'inside'));
% outsideTempFld  = find(contains(lower(tb.Properties.VariableNames), 'outside'));
inTemp          = mean(tb.(tb.Properties.VariableNames{insideTempFld}));
% outTemp         = mean(tb.(tb.Properties.VariableNames{outsideTempFld}));

title(sprintf('Economy (total drives:%1.0fkm)\nInside temperature:  %1.1fC\nOutside temperature: %1.1fC', sum(kmDeltaX), inTemp, outTemp))
legend show

end

if ismember('TimestampIST', tb.Properties.VariableNames)
    filename = sprintf('%s_%1.0fkm.png', datestr(tb.TimestampIST(1), 'yyyy_mm_dd_HHMM'), range(tb.Odometerkm));
elseif ismember('TimestampIDT', tb.Properties.VariableNames)
    filename = sprintf('%s_%1.0fkm.png', datestr(tb.TimestampIDT(1), 'yyyy_mm_dd_HHMM'), range(tb.Odometerkm));
end


if issave
    saveas(f, filename);
end
if isclose
    close all
end


end


%% service functions

function [en_loss, en_loss_c] = calcEnLoss(specs, vel, const_pwr, outTemp, intemp)

if nargin<3
    const_pwr = 0; % in Watts
end
if nargin<4
    outTemp = 22; % in Watts
end
if nargin<5
    intemp = 22; % in Watts
end
% wh/km!
Ts    	= [10:5:40];
rhos    = [1.24210000000000;1.21880000000000;1.19570000000000;1.17280000000000;1.14970000000000;1.12630000000000;1.10240000000000];

% deltaTemp = intemp - outTemp; % positive is to heat up
% 
% if deltaTemp>0  % heating needed?
%     deltaE = max(deltaTemp-8, 0)*150;
% else % cooling needed
%     deltaE = -deltaTemp*150;
% end
% fprintf('assumed constant operations power: %1.1fkWh', const_pwr/1e3)
% if deltaTemp>0 
%     fprintf(' + heating %1.1fkWh', deltaE/1e3);
% else
%     fprintf(' + cooling %1.1fkWh', deltaE/1e3);
% end
% fprintf('(deltaT~%1.0fC)\n', deltaTemp);
% const_pwr = const_pwr + deltaE;

% fprintf('cooling/heating + opertion adding %1.1fkW (deltaT=%1.0fC)\n', const_pwr/1e3, deltaTemp)

% mass[kg]
% vel[km/h]
% rho = GP.rho; % 1.21; % kg/m^3 % change to temperature dependant?
rho = interp1(Ts, rhos, outTemp, 'linear', 'extrap');

% rho = 1.204;
drag_coefficient = specs.drag_coeff;
front_cross_section = specs.front_cross_section;
rolling_resistance = specs.rolling_res; % tires+road
mass = specs.kg_vehicle_mass;
% drag
f_drag = front_cross_section*drag_coefficient*0.5*rho*(vel/3.6).^2; % N
enloss_drag = f_drag*100000/1000/3600*10; % Wh/km

% roll
enloss_roll = mass*9.81*rolling_resistance*100000/1000/3600*10; % Wh/km

en_loss = round(enloss_roll + enloss_drag);

ts = 1000./(vel/3.6); % 1 km using mps --> t[sec]
t = ts/60/60;
en_loss_c = const_pwr(:).*t(:);

en_loss = en_loss(:) + en_loss_c(:);
end

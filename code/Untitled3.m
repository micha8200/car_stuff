% plotDrive('C:\Users\micha\Documents\Python\Tesla\f09022024_0710_t09022024_2240.csv')
[f, p]=uigetfile('C:\Users\micha\Documents\Python\Tesla\*.csv');

if f==0
    return
end
tb = readtable(fullfile(p, f)); 

[ok, ixState]       = ismember(tb.ShiftState, {  'R' 'P' 'D'});
tb.ShiftState       = ixState - 2;

%%

tb.sTime        = datenum(tb.TimestampUTC)*86400;
tb.sTime        = int32(tb.sTime - tb.sTime(1));
tb.sDeltaTime   = 10*ones(size(tb.sTime));
tb.whDeltaE     = round(tb.PowerkW.*1e3.*tb.sDeltaTime * GP.wh_per_joule);

kwhUsed_trapz       = trapz(double(tb.sTime), tb.PowerkW)*GP.kwh_per_joule *1e3;
% kwhUsed_rect = sum(tb.PowerkW(1:end-1).*1e3.*diff(tb.sTime) * GP.wh_per_joule)/1e3;
kwhUsed_rect        = sum(tb.PowerkW.*1e3.*tb.sDeltaTime * GP.wh_per_joule)/1e3;
kwhUsed_batRng      = range(tb.BatteryRangekm) * GP.epa_wh_per_km/1e3;
kwhUsed_BatPct      = range(tb.BatteryLevelprc)/100 * GP.kwh_in_new_battery;


% tb.Speedkm_per_h = tb.Speedkm_per_h .* tb.ShiftState;

% count only 'D' or 'R' states, and add 10sec padding around them
ok                  = ismember(tb.ShiftState, [-1 1]);
tb.mDeltaBySpeed    = round(tb.Speedkm_per_h/3.6.*tb.sDeltaTime);
tb.mDeltaOdo        = round([diff(tb.Odometerkm);0]*1e3);
tb_dr               = tb(ok, :);
tb_park             = tb(~ok, :);
ixjump              = find(diff(tb_dr.sTime)>10);




fprintf('drive: %s\n', f)
fprintf('mean economy by power integration  : %1.3f wh/km\n', sum(tb_dr.whDeltaE)./range(tb_dr.Odometerkm))
fprintf('total energy used in drive         : %1.3f kWh\n', sum(tb_dr.whDeltaE)./1e3)
fprintf('total energy used in idle          : %1.3f kWh\n', sum(tb_park.whDeltaE)./1e3)
fprintf('total energy used in drive(by bat) : %1.3f kWh\n', range(tb_dr.BatteryRangekm.*GP.epa_wh_per_km/1e3))


mDeltaAlts  = diff(smooth(tb_dr.Elevationm, 0.05));
mDeltaAlts  = diff(tb_dr.Elevationm);
mTotalUps   = sum(mDeltaAlts(mDeltaAlts>0));
mTotalDwns  = sum(mDeltaAlts(mDeltaAlts<0));

kWhEnergyForUps = mTotalUps * GP.kwh_per_m_per_kg * GP.kg_vehicle_mass;
kWhEnergyForDowns = - 0.93*mTotalDwns * GP.kwh_per_m_per_kg * GP.kg_vehicle_mass;
kWhEnergyForUps - kWhEnergyForDowns

kWhGivenUp = 52.3 * GP.epa_wh_per_km / 1e3;
kWhGivenDwn = 46.1 * GP.epa_wh_per_km / 1e3;
kWhGivenUp - kWhGivenDwn

%%

ok = tb_dr.mDeltaBySpeed>0;
x = tb_dr.Odometerkm(ok); x = x - x(1);
y = tb_dr.whDeltaE(ok)./(tb_dr.mDeltaBySpeed(ok)/1e3);
y_sm = smooth(y, 0.1);


figure()
grid on
grid minor
hold on
% ok = ~isinf(wh2km) && ~isnan(wh2km);
% smy = smooth(tb_dr.Odometerkm(ok), wh2km(ok), 0.1);
plot(x, y_sm);
xlabel('odo[km]')
ylabel('econ [wh/km]')
ylim([-200 500])
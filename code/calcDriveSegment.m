function tb1 = calcDriveSegment(tb)

mass = 1750+150;

% tb1 = tb(drivesIxs{1},:);
% tb1 = tb(:, {'ShiftState' 'Elevationm' 'Speedkmperh' 'lineIx' 'OutsudeTempFromComfC'});
tb1 = tb;

tb1.kmDrive         = [0; diff(tb.Odometerkm)];
tb1.kwhEnergyUsed   = [0; diff(tb.BatteryRangekm)]*GP.epa_wh_per_km/1e3;
tb1.mHeightChange   = [0; diff(tb.Elevationm)];

%% filter out 100% long drive due to BMC

% delta(KIN+POT)~=0 but Energy still the same => discart these points!!!
if any(tb.BatteryLevelprc==100)
    Ekin    = 0.5*mass*(tb.Speedkmperh/3.6).^2;
    Epot    = mass*GP.g*tb.Elevationm;
    dE      = [diff(Ekin+Epot);0];
    dE_bat  = [diff(tb.BatteryRangekm);0];
    batteryRangeMaxkm = max(tb.BatteryRangekm(tb.BatteryLevelprc==100));
    isFrzBat = dE~=0 & dE_bat==0 & tb.BatteryRangekm==batteryRangeMaxkm;

    tb1(isFrzBat,:) = [];
end

end
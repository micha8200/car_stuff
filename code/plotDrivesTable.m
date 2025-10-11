function plotDrivesTable(tbDrives)
%% Energy Use + EPA

f(2)= figure();
f(2).Name = 'driveEnergy';


hold on
box on
grid on
grid minor

mass            = 1750+150;
% mass            = 1750;

kmDriven        = cumsum(tbDrives.kmDrive);
mHeightCum      = cumsum(tbDrives.mHeightChange);
E_REAL          = cumsum(tbDrives.kwhEnergyUsed);

E_EPA           = -GP.epa_wh_per_km*kmDriven/1e3;
Vel_mps         = tbDrives.Speedkmperh/3.6;

% E_real_r         = tb.BatteryRangekm*GP.epa_wh_per_km/1e3;
% kwhEnergyLevel  = loesmooth(E_real_r);
% kwhE0           = E_real_r(1);
% E_epa           = kwhE0 - GP.epa_wh_per_km*kmDist/1e3; % including height
kwhDeltaEnByHeight = mass*9.81*mHeightCum * GP.kwh_per_joule;
kwhCurrentSpeedEnergy = 0.5*mass*Vel_mps.^2*GP.kwh_per_joule;

% current height corrected EPA drive Energy state
E_EPA_HC        = E_EPA  - kwhDeltaEnByHeight; % EPA, only with reduction due to current height (E_pot)
E_REAL_HC       = E_REAL + kwhDeltaEnByHeight;

% current height and current velocity corrected EPA Energy drive state
E_EPA_HC_VC     = E_EPA-kwhDeltaEnByHeight-kwhCurrentSpeedEnergy;
E_REAL_HC_VC     = E_REAL+kwhDeltaEnByHeight+kwhCurrentSpeedEnergy;
E_REAL_HC_VC_sm  = smooth(kmDriven, E_REAL_HC_VC, 0.2);

plot(kmDriven, E_EPA, ':k', 'LineWidth', 2, 'DisplayName', 'E_{EPA} flat')
plot(kmDriven, E_EPA_HC, ':', 'LineWidth', 2, 'DisplayName', 'E_{EPA} w Epot', 'visible', 'off')
% plot(kmDist, E_EPA_HC_VC, ':', 'LineWidth', 2, 'DisplayName', 'E_{EPA} w Epot+Ekin')
plot(kmDriven, E_REAL, '-', 'LineWidth', 2, 'DisplayName', 'E_{cons}')
% plot(kmDist, E_REAL_HC, '-', 'LineWidth', 2, 'DisplayName', 'E_{battery} w Epot')
plot(kmDriven, E_REAL_HC_VC, '-', 'LineWidth', 2, 'DisplayName', 'E_{cons} K,P invariant', 'UserData', tbDrives, 'visible', 'off')
plot(kmDriven, E_REAL_HC_VC_sm, '-.', 'LineWidth', 2, 'DisplayName', 'E_{cons} K,P invariant smoothed', 'UserData', tbDrives)


title(sprintf('Energy state\ndifference from EPA:%1.2fkWh in %1.0fkm\nTotal efficiency: %1.0f Wh/km', E_EPA_HC_VC(end)-E_REAL(end), range(kmDriven), range(E_REAL*1e3)/range(kmDriven)))
ylabel('Energy [kWh]')
xlabel('Distace traveled [km]')
lg = legend('show');
lg.ItemHitFcn = @LegendItemHitFnc;

yyaxis right

plot(kmDriven, Vel_mps*3.6, '.-', 'displayname', 'velocity[km/h]')
legend show
ylabel('Velocity [km/h]')

legend show

xlim([0 max(kmDriven)])

end
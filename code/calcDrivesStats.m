function [econ_kwh_100km, econ_std_kwh_100km]=calcDrivesStats(tbDrives)
%% Energy Use + EPA


mass            = 1750+250;
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
E_EPA_HC        = E_EPA-kwhDeltaEnByHeight;
E_REAL_HC       = E_REAL + kwhDeltaEnByHeight;

% current height and current velocity corrected EPA Energy drive state
E_EPA_HC_VC     = E_EPA-kwhDeltaEnByHeight-kwhCurrentSpeedEnergy;
E_REAL_HC_VC     = E_REAL+kwhDeltaEnByHeight+kwhCurrentSpeedEnergy;

% sprintf('Energy state\ndifference from EPA:%1.2fkWh in %1.0fkm\nTotal efficiency: %1.0f Wh/km', E_EPA_HC_VC(end)-E_REAL(end), range(kmDriven), range(E_REAL*1e3)/range(kmDriven))

%
Npoints = length(E_REAL_HC_VC);
[P, S] = polyfit(kmDriven,E_REAL_HC_VC,1); % E*x + E0
Rinv = inv(S.R);
Pcov = sqrt(diag((Rinv*Rinv')*S.normr^2/S.df));


if nargout==0
    fprintf('estimated consumption(V~%1.0f km/hr, N=%1.0f): %1.1f ? %1.1f kWh/100km\n', mean(tbDrives.Speedkmperh), Npoints, P(1)*100, Pcov(1)*100);
else
    econ_kwh_100km = P(1)*100;
    econ_std_kwh_100km = Pcov(1)*100;
end
end
% use tessie site->settings->export(choose month)

temp_drive  = readtable('fulldata_exapmles\24_03\driving_states.csv');
temp_bat    = readtable('fulldata_exapmles\24_03\battery_states.csv');
temp_charg  = readtable('fulldata_exapmles\24_03\charging_states.csv');
temp_clim   = readtable('fulldata_exapmles\24_03\climate_states.csv');

temp_drive.Properties.VariableNames(:) = strrep(temp_drive.Properties.VariableNames', '_', '');
temp_bat.Properties.VariableNames(:) = strrep(temp_bat.Properties.VariableNames', '_', '');
temp_charg.Properties.VariableNames(:) = strrep(temp_charg.Properties.VariableNames', '_', '');
temp_clim.Properties.VariableNames(:) = strrep(temp_clim.Properties.VariableNames', '_', '');

%%

% tic
% t1 = cellfun(@timestamp2time, temp_drive.TimestampUTC(1:1000));
% toc
% 
% tic
% t2 = timestamp2time_array(char(temp_drive.TimestampUTC(1:1000)));
% toc

temp_drive.t    = timestamp2time_array(temp_drive.TimestampUTC);
temp_bat.t      = timestamp2time_array(temp_bat.TimestampUTC);
temp_clim.t     = timestamp2time_array(temp_clim.TimestampUTC);
temp_charg.t     = timestamp2time_array(temp_charg.TimestampUTC);

bg = find(diff(temp_charg.t)>0);
temp_drive.BatteryRangemi = interp1(temp_charg.t(bg), temp_charg.BatteryRangemi(bg), temp_drive.t, 'linear');
% temp_drive.interpEnRemainkWh = interp1(temp_bat.t(bg), temp_bat.EnergyRemainingkWh(bg), temp_drive.t, 'linear', 'extrap');

% temp_bat.EnergyRemainingkWh = 
% temp_drive.interpEnRemainkWh = temp_drive.interpBatLevel * 60/100;
figure();
hold on
% plot(temp_bat.t, temp_bat.EnergyRemainingkWh, 'o')
plot(temp_charg.t, temp_charg.BatteryRangemi, 'o')
plot(temp_drive.t, temp_drive.BatteryRangemi, '.')


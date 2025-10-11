% tbl = readtable('MY_drives\XP7YGCEK5SB651296-2025-08-11.csv');
% writetable(tbl, 'MY_drives\XP7YGCEK5SB651296_2025_08_11_simp.csv');

tbl1 = readtable('MY_drives\XP7YGCEK5SB651296_2025_08_11_simp.csv');


dx  = diff(tbl1.Odometer_km_);
dy  = diff(tbl1.EnergyRemaining_kWh_);
ok = [dx~=0;true];
ix0 = find(~ok);


tbl2 = tbl1(ok,:);

plot(tbl1.Odometer_km_(ok), tbl1.EnergyRemaining_kWh_(ok));

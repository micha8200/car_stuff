classdef GPY
    properties (Constant)
        g = 9.81
        kg_vehicle_mass = 1997+70
        kwh_per_joule   = 2.77778e-7
        wh_per_joule    = 2.77778e-4
        nstd            = 2 % 97% of population
        
        kwh_per_liter = 8.9
        kwh_per_m_per_kg = 2.72500218e-06
        hp_in_kw = 1.34102
        kwh_in_new_battery = 78.4
        
        epa_range = 526
        epa_wh_per_km   = 149.0 % 78.4/526.3
        wltp_range = 586
        wltp_wh_per_km   = 133.8 % 78.4/586
        
        
        
        drag_coeff = 0.22
        front_cross_section = 2.5 % 3: 2.6 Y:3.11
        rolling_res = 0.01
    end
end

%%

%{
drag_coefficient = 0.23;
front_cross_section = 2.3;
rolling_resistance = 0.011; % tires+road

%}
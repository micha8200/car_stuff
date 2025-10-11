classdef GP
    properties (Constant)
        g = 9.81
        kg_vehicle_mass = 1765
        kwh_per_joule   = 2.77778e-7
        wh_per_joule    = 2.77778e-4
        nstd            = 2 % 97% of population
        epa_wh_per_km   = 137 % 60/437
        kwh_per_liter = 8.9
        kwh_per_m_per_kg = 2.72500218e-06
        hp_in_kw = 1.34102
        kwh_in_new_battery = 57.6
        rho                   = 1.06; % 1.2 low temp  1.15 @35C/80%humidity   1.05 @35C/90%humidity;        
    end
end
function sepcs = tesla_models_specs(vin)
modelType       = vin(4);
modelYear       = strfind('MNPRSTVWXY', vin(10)) + 2020;
mc              = containers.Map({'S' 'K'}, {'RWD' 'AWD'});
motorConf       = mc(vin(8)); % S-single motor RWD K-sual motor AWD
sepcs.name   	= sprintf('Model%s-%1.0f %s', modelType, modelYear, motorConf);
def             = false;

if strcmpi(modelType,'3') && modelYear<=2023 && modelYear>=2021 && strcmpi(motorConf,'rwd')
    def = true;
    sepcs.kg_vehicle_mass       = 1765+70;
    sepcs.kwh_in_new_battery    = 60; % 60.5
    sepcs.drag_coeff            = 0.23;
    sepcs.front_cross_section   = 2.22; % 3: 2.6 Y:3.11;
    sepcs.rolling_res           = 0.0091; % Aero/Sport 0.0091  performance 0.011

    sepcs.epa_wh_per_km         = 137.8; % 60/437
    sepcs.epa_range             = 439; % 
    sepcs.wltp_wh_per_km      	= 123.2; % 60/491;
    sepcs.wltp_range          	= 491; 
elseif strcmpi(modelType,'y') && modelYear>=2025 && strcmpi(motorConf,'awd')
    def = true;
    sepcs.kg_vehicle_mass       = 1997+100;
    sepcs.kwh_in_new_battery    = 78.4;
    sepcs.drag_coeff            = 0.23; % 0.22
    sepcs.front_cross_section   = 2.65; % 2.6-2.7
    sepcs.rolling_res           = 0.009; %0.008 - 0.011

    sepcs.epa_wh_per_km         = 153.0; % 78.4/526.3 = 149.0
    sepcs.epa_range             = 526; % 526
    sepcs.wltp_wh_per_km      	= 133.8; % 78.4/586;
    sepcs.wltp_range          	= 586; % 586;
elseif strcmpi(modelType,'y') && modelYear>=2025 && strcmpi(motorConf,'rwd')
    def = true;
    sepcs.kg_vehicle_mass       = 1928+100;
    sepcs.kwh_in_new_battery    = 62.5;
    sepcs.drag_coeff            = 0.22; %0.22-0.23
    sepcs.front_cross_section   = 2.65; % 2.6-2.7
    sepcs.rolling_res           = 0.009; %0.008 - 0.011

    sepcs.epa_wh_per_km         = 143; % 62.5/440
    sepcs.epa_range             = 440; % 
    sepcs.wltp_wh_per_km      	= 125; % 62.5/500
    sepcs.wltp_range          	= 500; % 500
end
assert(def, sprintf('%1.0f model %s not found in specs', modelYear, modelType));
end
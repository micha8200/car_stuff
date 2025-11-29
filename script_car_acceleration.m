% string to measure car power

folder = 'E:\MyCodes\car_stuff\acc3';
folder = 'E:\MyCodes\car_stuff\accY';
% folder = 'E:\MyCodes\car_stuff\accTucson';
files = dir(fullfile(folder,'measure*.dragger'));
files = fullfile(folder, {files.name})';

% files = [
%        {'E:\MyCodes\car_stuff\MY_drives\accmeter\measure_202520101710.dragger'}
%     {'E:\MyCodes\car_stuff\MY_drives\accmeter\measure_202525101849.dragger'}
%     % {'E:\MyCodes\car_stuff\MY_drives\accmeter\measure_202529101800.dragger'}  
% ];


%  f=ma
% P = d/dt (E) = m*v*a
% a = P/m *1/v

mt = zeros(1, length(files)); % mean times 0-100
w2hp = 1.341 / 1e3;

if contains(folder, '3')
    P = 283 / w2hp; %1.341 * 1e3; % HP->kW -> W
    g_max_friction = 0.5; % 0.4-0.5 for RWD/FWD 0.8-0.9 for AWD
    m = 1714 + 100;
elseif contains(folder, 'Y')
    P = 300 / w2hp; %1.341 * 1e3; % HP->kW -> W
    g_max_friction = 0.85;   
    m =  1992 +100;
    
end
% figure(); t=0:0.1:10; plot(t, v_t_car_acc(t, P/m, g_max_friction)*3.6); grid on

vc = P/m/(9.81*g_max_friction);
tc = P/m/(9.81*g_max_friction).^2;

% later find both params using:
% a(v) =:
%
% v < v_c  & v > 20km/h  :   P/m*1/v_c   - constant and is actually equal to balance*mu_coeff
%
% v >= v_c & v < 100+km/k :   P/m*1/v     - diminishing acceleration as velocity gets higher
% above v_crit, engine can use peak power without slip
% v1_crit = min(v_crit(v_crit>0));
% v = v0+a*dt
% F*v = P
% a = P/(m*v)
% dv/dt = P/mv
% dv*v = P/m dt
% 0.5v^2|v1->v2 = P/m t| t1->t2
% 0.5(v2^2-v1^2) = P/m dt
% or in general: 0.5(v^2-v0^2)=P/m(t-t0)
% v = sqrt(v0^2+3P/m*(t-t0))


% figure()

cc = lines(7);
cdata = cell(length(files), 1);
for i=1:length(files) %1:length(files) 1 4 
    str  = file2string(files{i});
    stlines = strsplit(str, newline)';
    
    dataGPS     = readBlock(stlines, 'gpsData');
    if max(dataGPS.spd)<28 % (101km/h)
        continue
    end

    dataAcc     = readBlock(stlines, 'accRawData');
    dataAcc.a   = sqrt(dataAcc.x.^2+dataAcc.y.^2+dataAcc.z.^2);
    dataAcc.sT  = dataAcc.t/1e3;
    dataGPS.sT  = dataGPS.t/1e3;
    [~, ixmx]   = max(dataGPS.spd);
    ixmn        = find(dataGPS.spd<1, 1, 'last');
    dataGPS     = dataGPS(ixmn:ixmx,:); % dataGPS.spd(ixmn:ixmx)
    
    % extend first non-zero rise in velocity to the first sample, to allow
    % later to interpolate exact launch time
    dataGPS.spd(1)  = extrap_interp_start_t(dataGPS.sT, dataGPS.spd);
    if isinf(dataGPS.spd(1))
        continue
    end
    t0              = dataGPS.sT(1);
    te              = dataGPS.sT(end);
    dataGPS.sT = dataGPS.sT - t0;
    dataAcc.sT = dataAcc.sT - t0;

    % plot(dataGPS.sT, dataGPS.spd, '.')
    %%
    
    % find maximum of V(x)
    % find last 0 of V(x)
    
    t       = dataGPS.sT;
    v       = dataGPS.spd;
    tq          = (0:0.5:te-t0+1)';
    vq         = interp1(t, v, tq, 'spline');
    p           = fliplr(polyfit(tq, vq, 8)); % in m/s^2
    
    
    vq          = polyspan(tq, p);
    aq          = polyspan(tq, p, 1);
    
    t100        = pass_vel(tq, vq, 100/3.6);
    t0          = pass_vel(tq, vq, 0);
    mt(i)       = t100-t0;
    
    v_kmh           = dataGPS.spd*3.6;
    vq_kmh          = vq*3.6;
    t_s             = dataGPS.sT - t0;


    tq_s            = tq - t0;
    aq_g        = aq/9.81;
    t_max = max(tq);

    data1 = table();
    data1.ixfile = i*zeros(size(t_s));
    data1.t = t_s;
    data1.v_kmh = v_kmh;
    data1 = data1(data1.t>0, :);
    cdata{i} = data1;
    %%
    
    ic = mod(i-1, 7)+1;
    % ploy_t_v_a_curves(t_s, v_kmh/3.6, tq, vq, aq, cc(ic,:), P/m, g_max_friction);

end
%%

% v(t) = a*t at v<vc
% vc: 
data1 = vertcat(cdata{:});
data1 = sortrows(data1, 't');
p           = fliplr(polyfit(data1.t, data1.v_kmh/3.6, 10));
tq          = (0:0.1:max(data1.t))';
vq          = polyspan(tq, p);
aq          = polyspan(tq, p, 1);

ploy_t_v_a_curves(data1.t, data1.v_kmh/3.6, tq, vq, aq, cc(ic,:), P/m, g_max_friction);

%%

figure()
hold on
grid on

data2 = cdata{6};
plot(data2.t, data2.v_kmh, '.')
t=0:0.1:10; 

[v, p2m, u] = v_t_car_acc(t, P/m*0.97, 0.75);
t100 = calc_0_to_v(100/3.6, p2m, u);
plot([0 t100 t100],  [100 100 0], ':r', 'LineWidth',2)
plot(t, v*3.6); 
title(sprintf('pwr=%1.0fHP mu=%1.2f t_{0-100}=%1.2fsec', p2m*m*w2hp, u, t100))

% subplot(2, 2, 1); hold on
% plot(data1.t, data1.v_kmh, '.')
% plot(tq, vq*3.6)
% subplot(2, 2, 3)
% plot(tq, aq)
% subplot(2, 2, [2 4])
% hold on
% plot(vq*3.6, acc_func(vq, P, m, g_max_friction)/9.81, ':');
% plot(vq*3.6, aq/9.81)



%%
% subplot(2, 2, 1);
% title(sprintf('mean 0-100km/h: %1.2f sec', mean(mt)))

function tp = pass_vel(tq, vq, vtest)
ix1     = find(vq>vtest, 1, 'first');
tp      = interp1(vq(ix1-1:ix1), tq(ix1-1:ix1), vtest);
end

function vq = polyspan(tq, p, order)
% p: [p0 p1 p2...]
% tq: [N*1] query times
if nargin<3
    order = 0;
end
Np      = length(p);
for i=1:order
    p       = p(2:end).*(1:Np-1);
    Np      = length(p);
end
vddm    = tq.^(0:Np-1);
vq      = p*vddm';
end

function data = readBlock(lines, blockName)
lines1 = lines(contains(lines, blockName));
data = cell(length(lines1),1);
for i=1:length(lines1)
    st = strfind(lines1{i}, '{');
    en = strfind(lines1{i}, '}');
    data{i} = jsondecode(lines1{i}(st:en));
end
data = struct2table(vertcat(data{:}));
end
folder = 'E:\MyCodes\car_stuff\acc3';
% folder = 'E:\MyCodes\car_stuff\accY';
files = dir(fullfile(folder,'measure*.dragger'));
files = fullfile(folder, {files.name})';

w2hp = 1.341 / 1e3;

if contains(folder, '3')
    P = 280 / w2hp; %1.341 * 1e3; % HP->kW -> W
    g_max_friction = 0.5; % 0.4-0.5 for RWD/FWD 0.8-0.9 for AWD
    m = 1714 + 100;
    modeln = 'model3';
elseif contains(folder, 'Y')
    P = 320 / w2hp; %1.341 * 1e3; % HP->kW -> W
    g_max_friction = 0.85;   
    m =  1992 +100;
    modeln = 'modelY';
end
ix = 11;


str  = file2string(files{ix});
stlines = strsplit(str, newline)';

dataGPS     = readBlock(stlines, 'gpsData');
assert(max(dataGPS.spd)>28, 'at least 101 km/h must be reached!') % (101km/h)

dataAcc     = readBlock(stlines, 'accRawData');
dataAcc.a   = sqrt(dataAcc.x.^2+dataAcc.y.^2+dataAcc.z.^2);
dataAcc.sT  = dataAcc.t/1e3;
dataGPS.sT  = dataGPS.t/1e3;
dataGPS     = movevars(dataGPS, "sT", "Before", "spd");
[~, ixmx]   = max(dataGPS.spd);
ixmn        = find(dataGPS.spd<1, 1, 'last');
dataGPS     = dataGPS(ixmn:ixmx,:); % dataGPS.spd(ixmn:ixmx)

% extend first non-zero rise in velocity to the first sample, to allow
% later to interpolate exact launch time
dataGPS.spd(1)  = extrap_interp_start_t(dataGPS.sT, dataGPS.spd);
assert(~isinf(dataGPS.spd(1)))

te              = dataGPS.sT(end);
t0              = pass_vel(dataGPS.sT, dataGPS.spd, 0);
dataGPS.sT(1)   = t0;
dataGPS.spd(1)  = 0;
dataGPS.sT      = dataGPS.sT - t0;


figure()
hold on
grid on


plot(dataGPS.sT, dataGPS.spd*3.6, '.')


t=0:0.1:10; 

[v, p2m, u, tc]     = v_t_car_acc(t, P/m*0.85, 0.56);

ok      = dataGPS.spd<105/3.6;
vcce    = sqrt(sum((v_t_car_acc(dataGPS.sT(ok), p2m, u) - dataGPS.spd(ok)).^2)/(sum(ok)-1));

t100            = calc_0_to_v(100/3.6, p2m, u);
plot([0 t100 t100],  [100 100 0], ':r', 'LineWidth',2)
plot(t, v*3.6); 
yl = get(gca, 'ylim');
patch([0 tc tc 0],[0 0 1 1]*yl(2),[0 0 0], 'FaceAlpha', 0.1, 'edgecolor', 'none');
patch([tc t100 t100 tc],[0 0 1 1]*yl(2),[1 0 0], 'FaceAlpha', 0.1, 'edgecolor', 'none');
text(0, yl(2), sprintf('a_{const}=%1.2fg', u), 'VerticalAlignment','top')
text(tc, yl(2), sprintf('P_{const}=%1.1fhp', p2m*m*w2hp), 'VerticalAlignment','top')
title(sprintf('%s\nt_{0-100}=%1.2fsec\nSRMSE=%1.3f', modeln, t100, vcce))



%%

function tp = pass_vel(tq, vq, vtest)
ix1     = find(vq>vtest, 1, 'first');
tp      = interp1(vq(ix1-1:ix1), tq(ix1-1:ix1), vtest);
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

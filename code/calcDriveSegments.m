[f, p]=uigetfile('f*.csv');

if f>0
    fullfilename = fullfile(p, f);
    tb = readtable(fullfilename); 
    tb.Properties.VariableNames(:) = strrep(tb.Properties.VariableNames', '_', '');
    [ok, ixState]       = ismember(tb.ShiftState, {  'R' 'P' 'D'});
    tb.ShiftState       = ixState - 2; % 1:drive 0:park -1:reverse
end
%%

ntb = table();
copyfields = {
    'ShiftState'
    'Odometerkm'
    'BatteryLevelprc'
    'BatteryRangekm'
    'Elevationm'
    'Speedkmperh'
    'OutsideTempC'
    'InsideTempC'
    'Climate'
    };

ntb = tb(:, copyfields);
ntb.lineIx = [1:length(ntb.ShiftState)]';
ntb.TempDifferenceC = ntb.OutsideTempC - ntb.InsideTempC;
ntb.OutsudeTempFromComfC = ntb.OutsideTempC - 22;

isdrive = ntb.ShiftState==1;
shiftChanges = diff([0;isdrive]);
ixs = find(shiftChanges>0);
ixe = find(shiftChanges<0);

% drivesIxs = cell(size(ixs));
tbseg = cell(size(ixs));
for i=1:length(ixs)
%     drivesIxs{i} = ntb.lineIx(ixs(i):ixe(i));
    tbseg{i} = calcDriveSegment(ntb(ixs(i):ixe(i),:));
end

tbDrives = vertcat(tbseg{:});


answ = questdlg('Calculations:', 'title', 'statistics', 'plot whole drive', 'plot whole drive');

%%

if strcmpi(answ, 'plot whole drive')
    plotDrivesTable(tbDrives);
end

%% speed bins

if strcmpi(answ, 'statistics')



% bin index      1   2   3   4   5   6   7   8   9
% speed_bins  = [0  40  50  60  70  80  90  100 110 inf];
speed_bins  = [0  20 40   70  80 85 90 95 100 105 110 115];
speed_bins  = [speed_bins max(tbDrives.Speedkmperh)];
ixbin       = discretize(tbDrives.Speedkmperh, speed_bins);
tbDrives.ixSpeedBin     = ixbin;

Nbins = length(speed_bins)-1;
econs = zeros(1, Nbins);
econs_err = zeros(1, Nbins);
for i=1:Nbins
    [econs(i), econs_err(i)] = calcDrivesStats(tbDrives(ixbin==i,:));
end

meanSpeed   = accumarray(tbDrives.ixSpeedBin, tbDrives.Speedkmperh, [Nbins 1], @mean);
figure();
errorbar(meanSpeed, -econs, econs_err);
grid on

%% temperature bins

% % bin index                 1    2   3   4  
% tempDiffBins        = [      -2   4        ];
% % tempDiffBins        = [      0       ];
% tempDiffBins        = [min(tbDrives.OutsudeTempFromComfC) tempDiffBins max(tbDrives.OutsudeTempFromComfC)];
% ixbin               = discretize(tbDrives.OutsudeTempFromComfC, tempDiffBins);
% tbDrives.ixTempBin     = ixbin;

%%

% plotDrivesTable(tbDrives(ixbin==1,:));


%%

plotDrivesTable(tbDrives);

%% cal bin map


% [c, ic, icc]= unique(tbDrives{:, {'ixTempBin', 'ixSpeedBin'}}, 'rows');
% 
% ixc         = accumarray(icc, [1:length(icc)]', [], @(x) {x});
% econ_map    = zeros(length(tempDiffBins)-1, length(speed_bins)-1);
% for i=1:length(ixc)
%     econ_map(c(i,1), c(i,2)) = calcDrivesStats(tbDrives(ixc{i},:));
% end
% 
% econ_vals       = {};
% econ_bins       = [0 5 8 10 12 14 16 30];
% econ_bins_grid  = discretize(-econ_map, [econ_bins inf]);
% 
% % [x, y]=meshgrid(tempDiffBins(2:end), speed_bins(1:end-1));
% 
% figure();
% % surf(x', y', econ_map, econ_map);
% 
% 
% im = imagesc(econ_bins_grid);
% 
% customColorBar(econ_bins, true);
% % colormap jet;
% grid on
% ax = gca;
% meanTemp    = accumarray(tbDrives.ixTempBin, tbDrives.OutsudeTempFromComfC, [], @mean);
% meanSpeed   = accumarray(tbDrives.ixSpeedBin, tbDrives.Speedkmperh, [], @mean);
% 
% 
% ax.XTickLabel   = arrayfun(@num2str, round(meanSpeed), 'un', 0);
% ax.YTick        = 1:3;
% ax.YTickLabel   = arrayfun(@num2str, round(meanTemp), 'un', 0);


end

%%






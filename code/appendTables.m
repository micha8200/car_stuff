function appendTables()
%%

[f, p]=uigetfile('*.csv', 'MultiSelect',  'on');

c = cell(size(f));
for i=1:length(c)
    c{i} = readtable(fullfile(p, f{i}));
end


tt = vertcat(c{:});

tt = sortrows(tt, 'Timestamp_IDT_');

tt.Properties.VariableNames = strrep(tt.Properties.VariableNames, '_', '');

parts = strsplit(f{1}, '-');
vin = parts{1};

% year(tt.Timestamp_IDT_(1))
% month(tt.Timestamp_IDT_(1))
% day(tt.Timestamp_IDT_(1))

newf = fullfile(p, sprintf('%s_merged.csv', vin));

writetable(tt, newf)


end
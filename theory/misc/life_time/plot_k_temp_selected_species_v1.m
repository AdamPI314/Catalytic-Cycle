%% global settings
markers = {'+' , 'o' , '*' , 'x' , 'square' , 'diamond' , 'v' , '^' , '>' , '<' , 'pentagram' , 'hexagram' , '.'};

%% Current file directory
file_dir = fullfile(fileparts(mfilename('fullpath')), '..', '..', '..', '..', '..', '..', 'SOHR_DATA');
pic_dir = fullfile(fileparts(mfilename('fullpath')));

%% import time
fn_time = fullfile(file_dir, 'output', 'time_dlsode_M.csv');
delimiter = '';
formatSpec = '%f%[^\n\r]';
%% Open the text file.
fileID = fopen(fn_time,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN, 'ReturnOnError', false);

%% Close the text file.
fclose(fileID);
time_vec = dataArray{:, 1};
%% Clear temporary variables
clearvars fn_time delimiter formatSpec fileID dataArray ans;

%% import drc
fn_drc = fullfile(file_dir, 'output', 'drc_dlsode_M.csv');
delimiter = ',';
% For more information, see the TEXTSCAN documentation.
formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
%% Open the text file.
fileID = fopen(fn_drc,'r');
%% Read columns of data according to format string.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN, 'ReturnOnError', false);
%% Close the text file.
fclose(fileID);
%% Create output variable
drc_mat = [dataArray{1:end-1}];
%% Clear temporary variables
clearvars fn_conc delimiter formatSpec fileID dataArray ans;

%% chattering group drc
%% Initialize variables.
fn_chattering_drc = fullfile(file_dir, 'output', 'chattering_group_drc_dlsode_M.csv');

delimiter = ',';
formatSpec = '%f%f%f%f%f%f%f%f%f%[^\n\r]';
%% Open the text file.
fileID = fopen(fn_chattering_drc,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN, 'ReturnOnError', false);
%% Close the text file.
fclose(fileID);
%% Create output variable
chattering_drc_mat = [dataArray{1:end-1}];
%% Clear temporary variables
clearvars fn_chattering_drc delimiter formatSpec fileID dataArray ans;

%% import temperature
fn_temp = fullfile(file_dir, 'output', 'temperature_dlsode_M.csv');
delimiter = '';
formatSpec = '%f%[^\n\r]';
%% Open the text file.
fileID = fopen(fn_temp,'r');
%% Read columns of data according to format string.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN, 'ReturnOnError', false);
%% Close the text file.
fclose(fileID);
%% Allocate imported array to column variable names
temp_vec = dataArray{:, 1};
%% Clear temporary variables
clearvars fn_temp delimiter formatSpec fileID dataArray ans;

%% plot
fig = figure();
% https://www.mathworks.com/help/matlab/graphics_transition/why-are-plot-lines-different-colors.html
% https://www.mathworks.com/help/matlab/creating_plots/customize-graph-with-two-y-axes.html
%%
co = [    
    1   0   0
    ]; % 9th plot, red
set(fig,'defaultAxesColorOrder',co);
tau = 0.777660157519;
end_t = 0.9;

% plot chattering group drc
spe_idx_vec = [79, 76, 17, 28, 51, 49, 81, 77, 59, 95, 69, 63, 44];
N = length(spe_idx_vec);
for i = 1:N
    spe_idx_vec(i) = spe_idx_vec(i) + 1;
end
drc_mat = drc_mat(:, spe_idx_vec);

% sort by the reaction rates around 0.5 tau, idx == 3550 for example
sort_axis = round(0.6 * length(time_vec));
[B,I] = sort(1./drc_mat(sort_axis, :),'descend');

% legend name
old_legend_str = cell(N,1);
old_legend_str{1, 1} = 'nROOH';
old_legend_str{2, 1} = 'nRO';
old_legend_str{3, 1} = 'CH_2O';
old_legend_str{4, 1} = 'CH_3OOH';
old_legend_str{5, 1} = 'CH_3CH_2OOH';
old_legend_str{6, 1} = 'ethoxy';
old_legend_str{7, 1} = 'iROOH';
old_legend_str{8, 1} = 'iRO';
old_legend_str{9, 1} = 'C_3H_6';
old_legend_str{10, 1} = 'prod_2';
old_legend_str{11, 1} = 'acrolein';
old_legend_str{12, 1} = 'allyl';
old_legend_str{13, 1} = 'acetaldehyde';

legend_str = cell(N,1);

colors = lines(N);

delta_n = 1000;
H = gobjects(N);
for idx=1:N
%     spe_idx = ch_idx_vec(idx) + 1;
%     spe_idx = spe_idx_vec(I(idx));
    spe_idx = I(idx);
    legend_str{idx, 1} = old_legend_str{I(idx), 1};

    H(idx) = semilogy(time_vec, 1./drc_mat(:, spe_idx),...
        'LineWidth', 2, 'color', colors(idx, :), ...
        'HandleVisibility','off');

    hold on;
%     set(gca,'yscale','log');
    scatter(time_vec(1:delta_n:end), 1./drc_mat(1:delta_n:end, spe_idx),...
        'LineWidth', 2, 'MarkerEdgeColor', colors(idx, :), 'marker', markers{1, mod(idx-1, length(markers))+ 1}, ...
        'HandleVisibility','off');
    hold on;
    plot(nan, nan, 'LineWidth', 2, 'color', colors(idx, :), 'marker', markers{1, mod(idx-1, length(markers))+ 1});
    hold on;
end

%% conc
set(gca,'GridLineStyle','--');
xlabel('Time (seconds)', 'FontSize', 20);
ylabel('k^{-1} (seconds)', 'FontSize', 20);
ylim([10^-10, 10^5.35]);
% ylim([10^-12, 10^8]);

%% temperature
yyaxis right;
delta_n2 = 1000;
plot(time_vec, temp_vec, 'LineWidth', 2, 'color', 'r'); hold on;
pt = scatter(time_vec(1:delta_n2:end), temp_vec(1:delta_n2:end), 'MarkerEdgeColor', 'r', 'HandleVisibility','off');
ylabel('T (K)', 'FontSize', 20);
% set(gca, 'ytick', []);

%% global settings
grid on;
xlim([0, tau*end_t]);
% leg_h = legend(old_legend_str);
leg_h = legend(legend_str);
set(leg_h, 'FontSize', 11, 'Box', 'off');
set(leg_h, 'Location', 'North');



%% save to file
figname = strcat('selected_species_v1', '.png');
print(fig, fullfile(pic_dir, figname), '-r200', '-dpng');



%% global settings
fig_prefix = 'lifetime_low_T_500K_cv';
tau = 1899.7449542644167;
begin_t = 0.1;
end_t = 0.99;

markers = {'+' , 'o' , '*' , 'x' , 'square' , 'diamond' , 'v' , '^' , '>' , '<' , 'pentagram' , 'hexagram' , '.'};

%% Current file directory
file_dir = fullfile(fileparts(mfilename('fullpath')), '..', '..', '..', '..', '..', '..', '..', 'SOHR_DATA');
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

% plot chattering group drc
ch_idx_vec = [0, 1, 2, 3, 4, 5, 6, 7, 8];
N = length(ch_idx_vec);
colors = lines(N);

% sort by the reaction rates around 0.5 tau, idx == 3550 for example
sort_axis = round(0.1 * length(time_vec));
[B,I] = sort(1./chattering_drc_mat(sort_axis, :),'descend');

% legend name
old_legend_str = cell(N,1);
old_legend_str{1, 1} = 'CH_3 & CH_3OO';
old_legend_str{2, 1} = 'C_2H_5 & C_2H_5OO';
old_legend_str{3, 1} = 'acetyl & acetylperoxy';
old_legend_str{4, 1} = 'nR & nROO & QOOH_1 &O_2QOOH_1';
old_legend_str{5, 1} = 'iR & iROO';
old_legend_str{6, 1} = 'allyloxy & vinoxylmethyl';
old_legend_str{7, 1} = 'C_2H_4OH & O_2C_2H_4OH';
old_legend_str{8, 1} = 'QOOH_2 & O_2QOOH_2';
old_legend_str{9, 1} = 'QOOH_3 & O_2QOOH_3';

legend_str = cell(N,1);

delta_n = 8500;
H = gobjects(N);
for idx=1:N
%     group_idx = ch_idx_vec(idx) + 1;
    group_idx = ch_idx_vec(I(idx)) + 1;
    legend_str{idx, 1} = old_legend_str{I(idx), 1};
    if group_idx == 3
        init_idx = 350;
        H(idx) = semilogy(time_vec(init_idx:end), 1./chattering_drc_mat(init_idx:end, group_idx),...
        'LineWidth', 2, 'color', colors(idx, :), ...
        'HandleVisibility','off'); 
    else
        H(idx) = semilogy(time_vec, 1./chattering_drc_mat(:, group_idx),...
        'LineWidth', 2, 'color', colors(idx, :), ...
        'HandleVisibility','off'); 
    end
    hold on;
%     set(gca,'yscale','log');
    scatter(time_vec(1:delta_n:end), 1./chattering_drc_mat(1:delta_n:end, group_idx),...
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
ylim([10^-2.1, 10^0.1]);

%% temperature
yyaxis right;
plot(time_vec, temp_vec, 'LineWidth', 2, 'color', 'r'); hold on;
pt = scatter(time_vec(1:delta_n:end), temp_vec(1:delta_n:end), 'MarkerEdgeColor', 'r', 'HandleVisibility','off');
ylabel('T (K)', 'FontSize', 20);
% set(gca, 'ytick', []);

%% global settings
grid on;
xlim([tau*begin_t, tau*end_t]);
leg_h = legend(legend_str);
set(leg_h, 'FontSize', 11, 'Box', 'off');
% set(leg_h, 'Location', 'West');


%% save to file
figname = strcat(fig_prefix, '_all_chattering_group', '.png');
print(fig, fullfile(pic_dir, figname), '-r200', '-dpng');



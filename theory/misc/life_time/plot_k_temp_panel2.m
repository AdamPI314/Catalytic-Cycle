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


%% plot
fig = figure();
% https://www.mathworks.com/help/matlab/graphics_transition/why-are-plot-lines-different-colors.html
% https://www.mathworks.com/help/matlab/creating_plots/customize-graph-with-two-y-axes.html
co = [    0    0.4470    0.7410 % 1th plot
    0.8500    0.3250    0.0980 % 2nd plot
    0.9290    0.6940    0.1250 % 3rd plot
    0.4940    0.1840    0.5560 % 4th plot
    0.4660    0.6740    0.1880 % 5th plot
%     0.3010    0.7450    0.9330 % 6th plot
%     0.6350    0.0780    0.1840 % 7th plot
%     0         0    1.0000 % 8th plot, blue
    1   0   0]; % 9th plot, red
set(fig,'defaultAxesColorOrder',co)
%%
tau = 0.777660157519;
end_t = 0.9;

% plot drc
% npropyl
p1 = semilogy(time_vec, 1./drc_mat(:, 61), 'LineWidth', 2); hold on;
% npropyloo
p2 = semilogy(time_vec, 1./drc_mat(:, 79), 'LineWidth', 2); hold on;
% QOOH_1
p3 = semilogy(time_vec, 1./drc_mat(:, 88), 'LineWidth', 2); hold on;
% well_1
p4 = semilogy(time_vec, 1./drc_mat(:, 91), 'LineWidth', 2); hold on;
% chattering drc
ch_idx = 4;
p5 = semilogy(time_vec, 1./chattering_drc_mat(:, ch_idx), 'LineWidth', 2); hold on;

%% conc
set(gca,'GridLineStyle','--');
xlabel('Time (seconds)', 'FontSize', 20);
ylabel('Lifetime (seconds)', 'FontSize', 20);
ylim([10^-8.5, 10^-2.8]);

% %% temp
% yyaxis right
% delta_n = 1000;
% plot(time_vec, temp_vec, 'LineWidth', 2, 'color', 'r'); hold on;
% pt = scatter(time_vec(1:delta_n:end), temp_vec(1:delta_n:end), 'MarkerEdgeColor', 'r');
% ylabel('T (K)', 'FontSize', 20);
% % set(gca, 'ytick', []);

%% global settings
grid on;
xlim([0, tau*end_t]);
leg_h = legend([p1; p2; p3; p4; p5],'nR','nROO','QOOH_1','O_2QOOH_1','Chattering Group');
set(leg_h, 'FontSize', 16, 'Box', 'off');
set(leg_h, 'Location', 'North')


%% save to file
figname = strcat('drc_panel2', '.png');
print(fig, fullfile(pic_dir, figname), '-r200', '-dpng');



%% Current file directory
file_dir = fullfile(fileparts(mfilename('fullpath')), '..', '..', '..', '..', '..', '..', 'SOHR_DATA');

%% dlsode time
filename = fullfile(file_dir, 'output', 'time_dlsode_fraction.csv');

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

%% import concentration
fn_conc = fullfile(file_dir, 'output', 'concentration_dlsode_M.csv');
delimiter = ',';
% For more information, see the TEXTSCAN documentation.
formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
%% Open the text file.
fileID = fopen(fn_conc,'r');
%% Read columns of data according to format string.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN, 'ReturnOnError', false);
%% Close the text file.
fclose(fileID);
%% Create output variable
conc_mat = [dataArray{1:end-1}];
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

%% plot
fig = figure();
% https://www.mathworks.com/help/matlab/graphics_transition/why-are-plot-lines-different-colors.html
% https://www.mathworks.com/help/matlab/creating_plots/customize-graph-with-two-y-axes.html
co = [    0    0.4470    0.7410 % 1th plot
    0.8500    0.3250    0.0980 % 2nd plot
    0.9290    0.6940    0.1250 % 3rd plot
    0.4940    0.1840    0.5560 % 4th plot
    0.4660    0.6740    0.1880 % 5th plot
    0.3010    0.7450    0.9330 % 6th plot
    0.6350    0.0780    0.1840 % 7th plot
    0         0    1.0000 % 8th plot, blue
    1   0   0]; % 9th plot, red
set(fig,'defaultAxesColorOrder',co)
%%
tau = 0.777660157519;
end_t = 0.9;

% plot conc
% C3H6
p1 = semilogy(time_vec, conc_mat(:, 60), 'LineWidth', 2); hold on;
% CH2O
p2 = semilogy(time_vec, conc_mat(:, 18), 'LineWidth', 2); hold on;
% CO
p3 = semilogy(time_vec, conc_mat(:, 15), 'LineWidth', 2); hold on;
% prod_1
p4 = semilogy(time_vec, conc_mat(:, 95), 'LineWidth', 2); hold on;
% ipropyloo
p5 = semilogy(time_vec, conc_mat(:, 81), 'LineWidth', 2); hold on;
% acetaldehyde
p6 = semilogy(time_vec, conc_mat(:, 45), 'LineWidth', 2); hold on;
% propoxide
p7 = semilogy(time_vec, conc_mat(:, 87), 'LineWidth', 2); hold on;
% C2H4
p8 = semilogy(time_vec, conc_mat(:, 39), 'LineWidth', 2); hold on;


%% conc
set(gca,'GridLineStyle','--');
xlabel('Time (second)', 'FontSize', 20);
ylabel('[X] (mole\cdotL^{-1})', 'FontSize', 20);
ylim([10^-16, 10^-5.9]);

%% temp
yyaxis right
delta_n = 100;
plot(time_vec, temp_vec, 'LineWidth', 2, 'color', 'r'); hold on;
pt = scatter(time_vec(1:delta_n:end), temp_vec(1:delta_n:end), 'MarkerEdgeColor', 'r');
ylabel('T (K)', 'FontSize', 20);
% set(gca, 'ytick', []);

%% global settings
grid on;
xlim([0, tau*end_t]);
leg_h = legend([p1; p2; p3; p4; p5; p6; p7; p8],'C_3H_6','CH_2O','CO','OQ^{\prime}OOH_1','iROO','CH_3CHO','Propoxide', 'C_2H_4');
set(leg_h, 'FontSize', 14, 'Box', 'off');
set(leg_h, 'Location', 'NorthWest')


%% save to file
figname = strcat('panel1', '.png');
print(fig, fullfile(file_dir, 'output', figname), '-r200', '-dpng');



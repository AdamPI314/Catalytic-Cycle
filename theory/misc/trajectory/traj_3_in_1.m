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
    1   0   0]; % 8th plot
set(fig,'defaultAxesColorOrder',co)
%%
tau = 0.777660157519;
end_t = 0.258;

subplot(3, 1, 1);
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
% ipropylooh
p6 = semilogy(time_vec, conc_mat(:, 82), 'LineWidth', 2); hold on;
% acetaldehyde
p7 = semilogy(time_vec, conc_mat(:, 45), 'LineWidth', 2); hold on;

%% conc
set(gca,'GridLineStyle','--');
% xlabel('Time (second)', 'FontSize', 20);
ylabel('[X] (mole\cdotL^{-1})', 'FontSize', 20);
ylim([10^-16, 10^-9]);

%% temp
yyaxis right
delta_n = 100;
plot(time_vec, temp_vec, 'LineWidth', 2, 'color', 'r'); hold on;
scatter(time_vec(1:delta_n:end), temp_vec(1:delta_n:end), 'MarkerEdgeColor', 'r'); hold on;
ylabel('T (K)', 'FontSize', 20);
% set(gca, 'ytick', []);

%% global settings
grid on;
xlim([0, tau*end_t]);
leg_h = legend([p1; p2; p3; p4; p5; p6; p7],'C_3H_6','CH_2O','CO','OQ^{\prime}OOH','iROO','iROOH','CH_3CHO');
set(leg_h, 'FontSize', 14, 'Box', 'off');
set(leg_h, 'Location', 'NorthWest')

subplot(3, 1, 2);
% plot conc
% npropyl
p1 = semilogy(time_vec, conc_mat(:, 61), 'LineWidth', 2); hold on;
% npropyloo
p2 = semilogy(time_vec, conc_mat(:, 79), 'LineWidth', 2); hold on;
% QOOH_1
p3 = semilogy(time_vec, conc_mat(:, 88), 'LineWidth', 2); hold on;
% well_1
p4 = semilogy(time_vec, conc_mat(:, 91), 'LineWidth', 2); hold on;
% prod_1
p5 = semilogy(time_vec, conc_mat(:, 95), 'LineWidth', 2); hold on;
% frag_1
p6 = semilogy(time_vec, conc_mat(:, 102), 'LineWidth', 2); hold on;
% vinoxy
p7 = semilogy(time_vec, conc_mat(:, 47), 'LineWidth', 2); hold on;

%% conc
set(gca,'GridLineStyle','--');
% xlabel('Time (second)', 'FontSize', 20);
ylabel('[X] (mole\cdotL^{-1})', 'FontSize', 20);
ylim([10^-25, 10^-10]);

%% temp
yyaxis right
delta_n = 100;
plot(time_vec, temp_vec, 'LineWidth', 2, 'color', 'r'); hold on;
scatter(time_vec(1:delta_n:end), temp_vec(1:delta_n:end), 'MarkerEdgeColor', 'r'); hold on;
ylabel('T (K)', 'FontSize', 20);
% set(gca, 'ytick', []);

%% global settings
grid on;
xlim([0, tau*end_t]);
leg_h = legend([p1; p2; p3; p4; p5; p6; p7],'nR','nROO','QOOH_1','HOOQ^{\prime}OOH','OQ^{\prime}OOH','OQ^{\prime}O','CH_2CHO');
set(leg_h, 'FontSize', 14, 'Box', 'off');
set(leg_h, 'Location', 'NorthWest')

subplot(3,1,3);
% plot conc
% H2O
p1 = semilogy(time_vec, conc_mat(:, 12), 'LineWidth', 2); hold on;
% H2O2
p2 = semilogy(time_vec, conc_mat(:, 14), 'LineWidth', 2); hold on;
% HO2
p3 = semilogy(time_vec, conc_mat(:, 13), 'LineWidth', 2); hold on;
% OH
p4 = semilogy(time_vec, conc_mat(:, 11), 'LineWidth', 2); hold on;
% O
p5 = semilogy(time_vec, conc_mat(:, 9), 'LineWidth', 2); hold on;
% H
p6 = semilogy(time_vec, conc_mat(:, 4), 'LineWidth', 2); hold on;

%% conc
set(gca,'GridLineStyle','--');
xlabel('Time (second)', 'FontSize', 20);
ylabel('[X] (mole\cdotL^{-1})', 'FontSize', 20);
ylim([10^-25, 10^-9]);

%% temp
yyaxis right
delta_n = 100;
plot(time_vec, temp_vec, 'LineWidth', 2, 'color', 'r'); hold on;
scatter(time_vec(1:delta_n:end), temp_vec(1:delta_n:end), 'MarkerEdgeColor', 'r');
ylabel('T (K)', 'FontSize', 20);

%% global settings
grid on;
xlim([0, tau*end_t]);
leg_h = legend([p1; p2; p3; p4; p5; p6],'H_2O','H_2O_2','HO_2','OH','O','H');
set(leg_h, 'FontSize', 14, 'Box', 'off');
set(leg_h, 'Location', 'NorthWest')


%% save to file
figname = strcat('traj_3_in_1', '.png');
print(fig, fullfile(file_dir, 'output', figname), '-r200', '-dpng');
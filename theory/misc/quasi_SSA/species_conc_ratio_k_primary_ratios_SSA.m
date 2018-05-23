plot_exact_ssa_ratio = true;

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

%% import reaction rate
fn_R = fullfile(file_dir, 'output', 'reaction_rate_dlsode_M.csv');
delimiter = ',';
% For more information, see the TEXTSCAN documentation.
formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(fn_R,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN, 'ReturnOnError', false);
%% Close the text file.
fclose(fileID);
%% Create output variable
reaction_R_mat = [dataArray{1:end-1}];
%% Clear temporary variables
clearvars fn_R delimiter formatSpec fileID dataArray ans;

%% import SSA concentration
%% Initialize variables.
% fn_SSA = 'D:\Github\SOHR\projects\catalytic_cycle\theory\misc\quasi_SSA\data\chattering_group_ss_prob_dlsode_M.csv';
file_dir_SSA = fullfile(fileparts(mfilename('fullpath')));
fn_SSA = fullfile(file_dir_SSA, 'data', 'chattering_group_ss_prob_dlsode_M.csv');

delimiter = ',';
formatSpec = '%f%f%f%f%[^\n\r]';
%% Open the text file.
fileID = fopen(fn_SSA,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN, 'ReturnOnError', false);
%% Close the text file.
fclose(fileID);
SSA_vec = [dataArray{1:end-1}];
%% Clear temporary variables
clearvars fn_SSA delimiter formatSpec fileID dataArray ans;

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
%     0.6350    0.0780    0.1840 % 7th plot
%     0         0    1.0000 % 8th plot, blue
    1   0   0 % speceholder
    1   0   0 % speceholder
    1   0   0 % speceholder
    1   0   0 % speceholder
    1   0   0 % speceholder
    1   0   0 % speceholder
    1   0   0]; % 9th plot, red
set(fig,'defaultAxesColorOrder',co)
%%
tau = 0.777660157519;
end_t = 0.9;

SSA_delta = 250;
% plot conc
% npropyloo/npropyl
ratio1 = conc_mat(:, 79)./conc_mat(:, 61);
p1 = semilogy(time_vec, ratio1, 'color',co(1, :), 'LineWidth', 2); hold on;
ratio1_SSA= SSA_vec(:, 2)./SSA_vec(:, 1);
semilogy(time_vec, ratio1_SSA, ...
    'linestyle', ':', 'color',co(1, :), 'LineWidth', 2); hold on;
p1_SSA = scatter(time_vec(1:SSA_delta:end), ratio1_SSA(1:SSA_delta:end), ...
    'MarkerEdgeColor',co(1, :), 'Marker', 'x', 'LineWidth', 1.5); hold on;
% R_npropyloo * [npropyl] / (R_npropyl* [npropyloo])
p1_K = semilogy(time_vec, (reaction_R_mat(:, 1069+1).*conc_mat(:, 79)) ./ (reaction_R_mat(:, 1068+1).*conc_mat(:, 61)) , ...
    'color',co(1, :), 'linestyle', '--', 'LineWidth', 2); hold on;
% npropyloo/QOOH_1
ratio2 = conc_mat(:, 79)./conc_mat(:, 88);
p2 = semilogy(time_vec, ratio2, 'color',co(2, :), 'LineWidth', 2); hold on;
ratio2_SSA = SSA_vec(:, 2)./SSA_vec(:, 3);
semilogy(time_vec, ratio2_SSA, ...
    'linestyle', ':', 'color',co(2, :), 'LineWidth', 2); hold on;
p2_SSA = scatter(time_vec(1:SSA_delta:end), ratio2_SSA(1:SSA_delta:end), ...
    'MarkerEdgeColor',co(2, :), 'Marker', 'x', 'LineWidth', 1.5); hold on;
% prod_1
p2_K = semilogy(time_vec, (reaction_R_mat(:, 1080+1).*conc_mat(:, 79)) ./ (reaction_R_mat(:, 1081+1).*conc_mat(:, 88)) , ...
    'color',co(2, :), 'linestyle', '--', 'LineWidth', 2); hold on;
% well_1/QOOH_1
ratio3 = conc_mat(:, 91)./conc_mat(:, 88);
p3 = semilogy(time_vec, ratio3, 'color',co(3, :), 'LineWidth', 2); hold on;
ratio3_SSA = SSA_vec(:, 4)./SSA_vec(:, 3);
semilogy(time_vec, ratio3_SSA, ...
    'linestyle', ':', 'color',co(3, :), 'LineWidth', 2); hold on;
p3_SSA = scatter(time_vec(1:SSA_delta:end), ratio3_SSA(1:SSA_delta:end), ...
    'MarkerEdgeColor',co(3, :), 'Marker', 'x', 'LineWidth', 1.5); hold on;
% acetaldehyde
p3_K = semilogy(time_vec, (reaction_R_mat(:, 1116+1).*conc_mat(:, 91)) ./ (reaction_R_mat(:, 1117+1).*conc_mat(:, 88)) , ...
    'color',co(3, :), 'linestyle', '--', 'LineWidth', 2); hold on;

%% conc
set(gca,'GridLineStyle','--');
xlabel('Time (seconds)', 'FontSize', 20);
ylabel('Ratios', 'FontSize', 20);
ylim([10^1.7, 10^4.95]);

%% temp
yyaxis right
delta_n = 1000;
plot(time_vec, temp_vec, 'LineWidth', 2, 'color', 'r'); hold on;
pt = scatter(time_vec(1:delta_n:end), temp_vec(1:delta_n:end), 'MarkerEdgeColor', 'r');
ylabel('T (K)', 'FontSize', 20);
% set(gca, 'ytick', []);

%% global settings
grid on;
xlim([0, tau*end_t]);
leg_h = legend([p1; p1_K; p1_SSA; p2; p2_K; p2_SSA; p3; p3_K; p3_SSA;],... 
    '[nROO]/[nR]','K_{eq}','SSA', '[nROO]/[QOOH_1]','K_{eq}', 'SSA', '[O_2QOOH_1]/[QOOH_1]','K_{eq}', 'SSA');
set(leg_h, 'FontSize', 13, 'Box', 'off');
set(leg_h, 'Location', 'North')


%% save to file
figname = strcat('species_conc_ratio_k_primary_ratio_SSA_linear_algebra', '.png');
print(fig, fullfile(pic_dir, figname), '-r200', '-dpng');


%% A new figure, for the ratio of exact solution and SSA solution
if plot_exact_ssa_ratio == true
    fig = figure();
    % https://www.mathworks.com/help/matlab/graphics_transition/why-are-plot-lines-different-colors.html
    % https://www.mathworks.com/help/matlab/creating_plots/customize-graph-with-two-y-axes.html
    co = [    0    0.4470    0.7410 % 1th plot
        0.8500    0.3250    0.0980 % 2nd plot
        0.9290    0.6940    0.1250 % 3rd plot
        %     0.4940    0.1840    0.5560 % 4th plot
        %     0.4660    0.6740    0.1880 % 5th plot
        %     0.3010    0.7450    0.9330 % 6th plot
        %     0.6350    0.0780    0.1840 % 7th plot
        %     0         0    1.0000 % 8th plot, blue
        1   0   0]; % 9th plot, red
    set(fig,'defaultAxesColorOrder',co)
    %%
    tau = 0.777660157519;
    end_t = 0.9;
    
    SSA_delta = 250;
    % plot conc
    % npropyloo/npropyl
    ratio1 = conc_mat(:, 79)./conc_mat(:, 61);
    p1_r = semilogy(time_vec, ratio1./ratio1_SSA, 'color',co(1, :), 'LineWidth', 2); hold on;
    p2_r = semilogy(time_vec, ratio2./ratio2_SSA, 'color',co(2, :), 'LineWidth', 2); hold on;
    p3_r = semilogy(time_vec, ratio3./ratio3_SSA, 'color',co(3, :), 'LineWidth', 2); hold on;
    %% conc
    set(gca,'GridLineStyle','--');
    xlabel('Time (seconds)', 'FontSize', 20);
    ylabel('Ratios', 'FontSize', 20);
    ylim([0.98, 1.02]);
    
    %% temp
    yyaxis right
    delta_n = 1000;
    plot(time_vec, temp_vec, 'LineWidth', 2, 'color', 'r'); hold on;
    pt = scatter(time_vec(1:delta_n:end), temp_vec(1:delta_n:end), 'MarkerEdgeColor', 'r');
    ylabel('T (K)', 'FontSize', 20);
    % set(gca, 'ytick', []);
    
    %% global settings
    grid on;
    xlim([0, tau*end_t]);
    leg_h = legend([p1_r; p2_r; p3_r;],...
        '[nROO]/[nR]','[nROO]/[QOOH_1]','[O_2QOOH_1]/[QOOH_1]');
    set(leg_h, 'FontSize', 13, 'Box', 'off');
    set(leg_h, 'Location', 'North')
    
    %% save to file
    figname = strcat('ratio_exact_over_SSA', '.png');
    print(fig, fullfile(pic_dir, figname), '-r200', '-dpng');
end
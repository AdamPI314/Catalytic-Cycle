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

%% import time
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
    1   0   0]; % 9th plot, red
set(fig,'defaultAxesColorOrder',co)
%%
tau = 0.777660157519;
end_t = 0.9;

% plot conc
% npropyloo/npropyl
p1 = semilogy(time_vec, conc_mat(:, 79)./conc_mat(:, 61), 'color',co(1, :), 'LineWidth', 2); hold on;
% R_npropyloo * [npropyl] / (R_npropyl* [npropyloo])
idx_1_top = [1069, 706];
for idx=1:length(idx_1_top)
    r_idx = idx_1_top(idx) + 1;
    if idx == 1
        R_1_top = reaction_R_mat(:, r_idx);
    else
        R_1_top = (R_1_top + reaction_R_mat(:, r_idx));
    end
end

idx_1_bottom = [707, 1068];
for idx=1:length(idx_1_bottom)
    r_idx = idx_1_bottom(idx) + 1;
    if idx == 1
        R_1_bottom = reaction_R_mat(:, r_idx);
    else
        R_1_bottom = (R_1_bottom + reaction_R_mat(:, r_idx));
    end
end
% R_1_bottom = reaction_R_mat(:, 1068+1);
p1_eq = semilogy(time_vec, (R_1_top.*conc_mat(:, 79)) ./ (R_1_bottom.*conc_mat(:, 61)) , ...
    'color',co(1, :), 'linestyle', '-.', 'LineWidth', 2); hold on;
% npropyloo/QOOH_1
p2 = semilogy(time_vec, conc_mat(:, 79)./conc_mat(:, 88), 'color',co(2, :), 'LineWidth', 2); hold on;
% prod_1
idx_2_top = [1080];
for idx=1:length(idx_2_top)
    r_idx = idx_2_top(idx) + 1;
    if idx == 1
        R_2_top = reaction_R_mat(:, r_idx);
    else
        R_2_top = (R_2_top + reaction_R_mat(:, r_idx));
    end
end
idx_2_bottom = [1081];
for idx=1:length(idx_2_bottom)
    r_idx = idx_2_bottom(idx) + 1;
    if idx == 1
        R_2_bottom = reaction_R_mat(:, r_idx);
    else
        R_2_bottom = (R_2_bottom + reaction_R_mat(:, r_idx));
    end
end
p2_eq = semilogy(time_vec, (R_2_top.*conc_mat(:, 79)) ./ (R_2_bottom.*conc_mat(:, 88)) , ...
    'color',co(2, :), 'linestyle', '-.', 'LineWidth', 2); hold on;
% well_1/QOOH_1
p3 = semilogy(time_vec, conc_mat(:, 91)./conc_mat(:, 88), 'color',co(3, :), 'LineWidth', 2); hold on;
% acetaldehyde
idx_3_top = [1116];
for idx=1:length(idx_3_top)
    r_idx = idx_3_top(idx) + 1;
    if idx == 1
        R_3_top = reaction_R_mat(:, r_idx);
    else
        R_3_top = (R_3_top + reaction_R_mat(:, r_idx));
    end
end
idx_3_bottom = [1117];
for idx=1:length(idx_3_bottom)
    r_idx = idx_3_bottom(idx) + 1;
    if idx == 1
        R_3_bottom = reaction_R_mat(:, r_idx);
    else
        R_3_bottom = (R_3_bottom + reaction_R_mat(:, r_idx));
    end
end
p3_eq = semilogy(time_vec, (R_3_top.*conc_mat(:, 91)) ./ (R_3_bottom.*conc_mat(:, 88)) , ...
    'color',co(3, :), 'linestyle', '-.', 'LineWidth', 2); hold on;

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
leg_h = legend([p1; p1_eq; p2; p2_eq; p3; p3_eq],'[nROO]/[nR]','K_{eq}','[nROO]/[QOOH_1]','K_{eq}','[O_2QOOH_1]/[QOOH_1]','K_{eq}');
set(leg_h, 'FontSize', 14, 'Box', 'off');
set(leg_h, 'Location', 'South')


%% save to file
figname = strcat('species_conc_ratio_k_primary_ratio_v2', '.png');
print(fig, fullfile(pic_dir, figname), '-r200', '-dpng');
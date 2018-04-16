%% global settings
spe_name = 'HO2';

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
%     0.3010    0.7450    0.9330 % 6th plot
%     0.6350    0.0780    0.1840 % 7th plot
%     0   0   1 % placeholder
%     0   0.5   0 % placeholder
%     0   0.75   0.75 % placeholder
%     0.7500         0    0.7500 % placeholder
%     0.7500    0.7500         0 % placeholder
%     0.2500    0.2500    0.2500 % placeholder
    1   0   0 % bl
    ]; 
set(fig,'defaultAxesColorOrder',co)
%%
tau = 0.777660157519;
end_t = 0.9;

%% plot
R_idx = [24, 27, 29, 31, 33, 35, 40, 42, 44, 51, 56, 67, 92, 95, 99, 103, 111, 113, 125, 143, 166, 169, 179, 188, 191, 214, 219, 243, 314, 319, 339, 342, 357, 364, 374, 380, 393, 414, 419, 437, 438, 486, 501, 506, 538, 540, 555, 557, 559, 590, 604, 607, 616, 630, 633, 659, 670, 687, 692, 724, 726, 741, 743, 821, 823, 825, 834, 836, 838, 853, 889, 903, 907, 909, 923, 925, 998, 1009, 1031, 1033, 1035, 1074, 1082, 1088, 1092, 1100, 1106, 1110, 1115, 1120, 1136, 1138, 1140, 1156, 1158, 1160, 1164, 1178, 1180, 1182, 1192, 1194, 1196, 1204, 1206, 1208];
for idx=1:length(R_idx)
    R_idx(idx) = R_idx(idx) + 1;
end

R_name = {'H+O_2(+M) \rightarrow HO_2(+M)', 'H_2+O_2 \rightarrow HO_2+H', 'OH+OH \rightarrow HO_2+H', 'O_2+OH \rightarrow HO_2+O', 'H_2O+O_2 \rightarrow HO_2+OH', 'H_2O_2+O_2 \rightarrow HO_2+HO_2', 'H_2O_2+H \rightarrow HO_2+H_2', 'H_2O_2+O \rightarrow OH+HO_2', 'H_2O_2+OH \rightarrow HO_2+H_2O', 'CO_2+OH \rightarrow CO+HO_2', 'HCO+O_2 \rightarrow CO+HO_2', 'CO_2+OH+H \rightarrow HCO+HO_2', 'CH_2O+O_2 \rightarrow HCO+HO_2', 'HCO+H_2O_2 \rightarrow CH_2O+HO_2', 'OCH_2OOH \rightarrow CH_2O+HO_2', 'HOCH_2OOH+O_2 \rightarrow HOCH_2OO+HO_2', 'CH_3O+OH \rightarrow CH_3+HO_2', 'CH_4+O_2 \rightarrow CH_3+HO_2', 'CH_3+H_2O_2 \rightarrow CH_4+HO_2', 'CH_3OOH+O_2 \rightarrow CH_3OO+HO_2', 'CH_2OH+O_2 \rightarrow CH_2O+HO_2', 'CH_2O+H_2O_2 \rightarrow CH_2OH+HO_2', 'HOCH_2O+OH \rightarrow CH_2OH+HO_2', 'CH_3O+O_2 \rightarrow CH_2O+HO_2', 'CH_2O+H_2O_2 \rightarrow CH_3O+HO_2', 'CH_3OH+O_2 \rightarrow CH_2OH+HO_2', 'CH_2OH+H_2O_2 \rightarrow CH_3OH+HO_2', 'CH_2O+OH \rightarrow CH_2+HO_2', 'C_2H_6+O_2 \rightarrow C_2H_5+HO_2', 'C_2H_5+H_2O_2 \rightarrow C_2H_6+HO_2', 'ethoxy+OH \rightarrow C_2H_5+HO_2', 'ethoxy+O_2 \rightarrow acetaldehyde+HO_2', 'CH_3CH_2OOH+O_2 \rightarrow CH_3CH_2OO+HO_2', 'C_2H_5+O_2 \rightarrow C_2H_4+HO_2', 'CH_3CH_2OO \rightarrow C_2H_4+HO_2', 'CH_2CH_2OOH \rightarrow C_2H_4+HO_2', 'oxiranyl+H_2O_2 \rightarrow oxirane+HO_2', 'acetaldehyde+O_2 \rightarrow acetyl+HO_2', 'acetyl+H_2O_2 \rightarrow acetaldehyde+HO_2', 'CH_3CO_3H+O_2 \rightarrow acetylperoxy+HO_2', 'H_2O_2+acetylperoxy \rightarrow HO_2+CH_3CO_3H', 'C_2H_4+O_2 \rightarrow C_2H_3+HO_2', 'oxirane+OH \rightarrow C_2H_4+HO_2', 'C_2H_3+O_2 \rightarrow HO_2+C_2H_2', 'ethanol+O_2 \rightarrow CH_2CH_2OH+HO_2', 'ethanol+O_2 \rightarrow CH_3CHOH+HO_2', 'CH_2CH_2OH+H_2O_2 \rightarrow ethanol+HO_2', 'CH_3CHOH+H_2O_2 \rightarrow ethanol+HO_2', 'ethoxy+H_2O_2 \rightarrow ethanol+HO_2', 'CH_3CHOH+O_2 \rightarrow acetaldehyde+HO_2', 'acetone+O_2 \rightarrow propen2oxy+HO_2', 'propen2oxy+H_2O_2 \rightarrow acetone+HO_2', 'ipropyloxy+O_2 \rightarrow acetone+HO_2', 'acrolein+O_2 \rightarrow CH_2CHCO+HO_2', 'CH_2CHCO+H_2O_2 \rightarrow acrolein+HO_2', 'propionyl+H_2O_2 \rightarrow propanal+HO_2', 'propanal+O_2 \rightarrow propionyl+HO_2', 'CH_3OCH_2+H_2O_2 \rightarrow CH_3OCH_3+HO_2', 'CH_3OCH_3+O_2 \rightarrow CH_3OCH_2+HO_2', 'C_3H_8+O_2 \rightarrow ipropyl+HO_2', 'C_3H_8+O_2 \rightarrow npropyl+HO_2', 'ipropyl+H_2O_2 \rightarrow C_3H_8+HO_2', 'npropyl+H_2O_2 \rightarrow C_3H_8+HO_2', 'allyl+H_2O_2 \rightarrow C_3H_6+HO_2', 'propen1yl+H_2O_2 \rightarrow C_3H_6+HO_2', 'propen2yl+H_2O_2 \rightarrow C_3H_6+HO_2', 'C_3H_6+O_2 \rightarrow allyl+HO_2', 'C_3H_6+O_2 \rightarrow propen1yl+HO_2', 'C_3H_6+O_2 \rightarrow propen2yl+HO_2', 'propen1ol+OH \rightarrow C_3H_6+HO_2', 'C_2H_4+HCO+OH \rightarrow propen1yl+HO_2', 'CH_3+ketene+OH \rightarrow propen2yl+HO_2', 'npropyloxy+OH \rightarrow npropyl+HO_2', 'ipropyloxy+OH \rightarrow ipropyl+HO_2', 'npropylooh+O_2 \rightarrow npropyloo+HO_2', 'ipropylooh+O_2 \rightarrow ipropyloo+HO_2', 'allyloxy+O_2 \rightarrow acrolein+HO_2', 'CH_2O+C_2H_3+H_2O_2 \rightarrow propen1ol+HO_2', 'prod_2 \rightarrow allyl+HO_2', 'acrolein+H_2O \rightarrow allyl+HO_2', 'allyloxy+OH \rightarrow allyl+HO_2', 'O_2+npropyl \rightarrow HO_2+C_3H_6', 'npropyloo \rightarrow HO_2+C_3H_6', 'QOOH_2 \rightarrow HO_2+C_3H_6', 'QOOH_1 \rightarrow HO_2+C_3H_6', 'O_2+ipropyl \rightarrow HO_2+C_3H_6', 'ipropyloo \rightarrow HO_2+C_3H_6', 'QOOH_3 \rightarrow HO_2+C_3H_6', 'OH+propoxide \rightarrow HO_2+C_3H_6', 'O_2+QOOH_1 \rightarrow HO_2+prod_2', 'O_2+QOOH_2 \rightarrow HO_2+prod_2', 'O_2+QOOH_2 \rightarrow HO_2+prod_6', 'O_2+QOOH_2 \rightarrow HO_2+prod_7', 'O_2+QOOH_3 \rightarrow HO_2+prod_2', 'O_2+QOOH_3 \rightarrow HO_2+prod_6', 'O_2+QOOH_3 \rightarrow HO_2+prod_7', 'well_1 \rightarrow HO_2+prod_2', 'well_2 \rightarrow HO_2+prod_2', 'well_2 \rightarrow HO_2+prod_6', 'well_2 \rightarrow HO_2+prod_7', 'well_3 \rightarrow HO_2+prod_2', 'well_3 \rightarrow HO_2+prod_6', 'well_3 \rightarrow HO_2+prod_7', 'well_5 \rightarrow HO_2+prod_2', 'well_5 \rightarrow HO_2+prod_6', 'well_5 \rightarrow HO_2+prod_7'};

R_mat = reaction_R_mat(:, R_idx);

% sort by the reaction rates around 0.5 tau, idx == 3550 for example
sort_axis = round(0.5 * length(time_vec));
[B,I] = sort(R_mat(sort_axis, :),'descend');

for idx=1:length(I)
    r_idx = I(idx);
    if idx == 1
        R_total = R_mat(:, r_idx);
    else
        R_total = (R_total + R_mat(:, r_idx));
    end
end

N = 5;

% graph handler
H = gobjects(N);

for idx=1:N
    r_idx = I(idx);
    H(idx) = semilogy(time_vec, (R_mat(:, r_idx)) ./ (R_total), ...
        'LineWidth', 2); hold on;
end



%% settings
set(gca,'GridLineStyle','--');
xlabel('Time (seconds)', 'FontSize', 20);
ylabel('Fraction', 'FontSize', 20);
ylim([10^-2, 10^0]);

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
leg_h = legend([H(1);H(2);H(3);H(4);H(5)], ...
    R_name{1,I(1)},R_name{1,I(2)},R_name{1,I(3)},R_name{1,I(4)},R_name{1,I(5)});
set(leg_h, 'FontSize', 12, 'Box', 'off');
set(leg_h, 'Location', 'South')


%% save to file
figname = strcat(spe_name, '_source_reaction_ratio_top_5', '.png');
print(fig, fullfile(pic_dir, figname), '-r200', '-dpng');
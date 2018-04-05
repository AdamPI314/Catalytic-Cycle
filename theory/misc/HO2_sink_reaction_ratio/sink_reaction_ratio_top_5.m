%% global settings
spe_name = 'HO2';

%% Current file directory
file_dir = fullfile(fileparts(mfilename('fullpath')), '..', '..', '..', '..', '..', '..', 'SOHR_DATA');

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
R_idx = [922, 1111, 45, 1193, 1205, 50, 556, 418, 94, 218, 558, 215, 725, 168, 539, 34, 727, 343, 671, 888, 415, 124, 98, 381, 1101, 25, 436, 824, 30, 500, 906, 820, 190, 1197, 1183, 606, 57, 1181, 554, 439, 392, 617, 1137, 178, 315, 110, 1114, 26, 1107, 318, 102, 1207, 32, 693, 631, 852, 742, 365, 902, 658, 507, 1165, 1083, 189, 142, 375, 1032, 1161, 167, 338, 66, 591, 1093, 242, 1141, 837, 28, 1195, 1157, 686, 1089, 356, 112, 1075, 839, 740, 1030, 41, 1179, 1121, 487, 1159, 1139, 999, 632, 605, 1209, 1008, 541, 908, 924, 43, 93, 822, 1034, 835];
for idx=1:length(R_idx)
    R_idx(idx) = R_idx(idx) + 1;
end

R_name = {'npropyloo+HO_2 \rightarrow npropylooh+O_2', 'HO_2+C_3H_6 \rightarrow QOOH_3', 'HO_2+H_2O \rightarrow H_2O_2+OH', 'HO_2+prod_2 \rightarrow well_3', 'HO_2+prod_2 \rightarrow well_5', 'CO+HO_2 \rightarrow CO_2+OH', 'ethanol+HO_2 \rightarrow CH_3CHOH+H_2O_2', 'acetaldehyde+HO_2 \rightarrow acetyl+H_2O_2', 'CH_2O+HO_2 \rightarrow HCO+H_2O_2', 'CH_3OH+HO_2 \rightarrow CH_2OH+H_2O_2', 'ethanol+HO_2 \rightarrow ethoxy+H_2O_2', 'CH_2OH+HO_2 \rightarrow CH_3OH+O_2', 'ipropyl+HO_2 \rightarrow C_3H_8+O_2', 'CH_2OH+HO_2 \rightarrow CH_2O+H_2O_2', 'CH_2CH_2OH+HO_2 \rightarrow ethanol+O_2', 'HO_2+HO_2 \rightarrow H_2O_2+O_2', 'npropyl+HO_2 \rightarrow C_3H_8+O_2', 'acetaldehyde+HO_2 \rightarrow ethoxy+O_2', 'propionyl+HO_2 \rightarrow propanal+O_2', 'propen1yl+HO_2 \rightarrow C_2H_4+HCO+OH', 'acetyl+HO_2 \rightarrow acetaldehyde+O_2', 'CH_4+HO_2 \rightarrow CH_3+H_2O_2', 'CH_2O+HO_2 \rightarrow OCH_2OOH', 'C_2H_4+HO_2 \rightarrow CH_2CH_2OOH', 'HO_2+C_3H_6 \rightarrow O_2+ipropyl', 'HO_2(+M) \rightarrow H+O_2(+M)', 'acetylperoxy+HO_2 \rightarrow CH_3CO_3H+O_2', 'C_3H_6+HO_2 \rightarrow propen2yl+H_2O_2', 'HO_2+O \rightarrow O_2+OH', 'C_2H_4+HO_2 \rightarrow oxirane+OH', 'npropyl+HO_2 \rightarrow npropyloxy+OH', 'C_3H_6+HO_2 \rightarrow allyl+H_2O_2', 'CH_3O+HO_2 \rightarrow CH_2O+H_2O_2', 'HO_2+prod_7 \rightarrow well_3', 'HO_2+prod_7 \rightarrow well_2', 'acetone+HO_2 \rightarrow propen2oxy+H_2O_2', 'CO+HO_2 \rightarrow HCO+O_2', 'HO_2+prod_6 \rightarrow well_2', 'ethanol+HO_2 \rightarrow CH_2CH_2OH+H_2O_2', 'HO_2+CH_3CO_3H \rightarrow H_2O_2+acetylperoxy', 'oxirane+HO_2 \rightarrow oxiranyl+H_2O_2', 'acetone+HO_2 \rightarrow ipropyloxy+O_2', 'HO_2+prod_2 \rightarrow O_2+QOOH_2', 'CH_2OH+HO_2 \rightarrow HOCH_2O+OH', 'C_2H_5+HO_2 \rightarrow C_2H_6+O_2', 'CH_3+HO_2 \rightarrow CH_3O+OH', 'HO_2+C_3H_6 \rightarrow OH+propoxide', 'HO_2+H \rightarrow H_2+O_2', 'HO_2+C_3H_6 \rightarrow ipropyloo', 'C_2H_6+HO_2 \rightarrow C_2H_5+H_2O_2', 'HOCH_2OO+HO_2 \rightarrow HOCH_2OOH+O_2', 'HO_2+prod_6 \rightarrow well_5', 'HO_2+OH \rightarrow H_2O+O_2', 'CH_3OCH_2+HO_2 \rightarrow CH_3OCH_3+O_2', 'CH_2CHCO+HO_2 \rightarrow acrolein+O_2', 'C_3H_6+HO_2 \rightarrow propen1ol+OH', 'C_3H_8+HO_2 \rightarrow npropyl+H_2O_2', 'C_2H_4+HO_2 \rightarrow C_2H_5+O_2', 'propen2yl+HO_2 \rightarrow CH_3+ketene+OH', 'propanal+HO_2 \rightarrow propionyl+H_2O_2', 'HO_2+C_2H_2 \rightarrow C_2H_3+O_2', 'HO_2+prod_2 \rightarrow well_1', 'HO_2+C_3H_6 \rightarrow npropyloo', 'CH_2O+HO_2 \rightarrow CH_3O+O_2', 'CH_3OO+HO_2 \rightarrow CH_3OOH+O_2', 'C_2H_4+HO_2 \rightarrow CH_3CH_2OO', 'allyl+HO_2 \rightarrow acrolein+H_2O', 'HO_2+prod_7 \rightarrow O_2+QOOH_3', 'CH_2O+HO_2 \rightarrow CH_2OH+O_2', 'C_2H_5+HO_2 \rightarrow ethoxy+OH', 'HCO+HO_2 \rightarrow CO_2+OH+H', 'acetaldehyde+HO_2 \rightarrow CH_3CHOH+O_2', 'HO_2+C_3H_6 \rightarrow QOOH_1', 'CH_2+HO_2 \rightarrow CH_2O+OH', 'HO_2+prod_7 \rightarrow O_2+QOOH_2', 'propen1yl+HO_2 \rightarrow C_3H_6+O_2', 'HO_2+H \rightarrow OH+OH', 'HO_2+prod_6 \rightarrow well_3', 'HO_2+prod_2 \rightarrow O_2+QOOH_3', 'CH_3OCH_3+HO_2 \rightarrow CH_3OCH_2+H_2O_2', 'HO_2+C_3H_6 \rightarrow QOOH_2', 'CH_3CH_2OO+HO_2 \rightarrow CH_3CH_2OOH+O_2', 'CH_3+HO_2 \rightarrow CH_4+O_2', 'HO_2+C_3H_6 \rightarrow O_2+npropyl', 'propen2yl+HO_2 \rightarrow C_3H_6+O_2', 'C_3H_8+HO_2 \rightarrow ipropyl+H_2O_2', 'allyl+HO_2 \rightarrow prod_2', 'HO_2+H_2 \rightarrow H_2O_2+H', 'HO_2+prod_2 \rightarrow well_2', 'HO_2+prod_2 \rightarrow O_2+QOOH_1', 'C_2H_3+HO_2 \rightarrow C_2H_4+O_2', 'HO_2+prod_6 \rightarrow O_2+QOOH_3', 'HO_2+prod_6 \rightarrow O_2+QOOH_2', 'acrolein+HO_2 \rightarrow allyloxy+O_2', 'acrolein+HO_2 \rightarrow CH_2CHCO+H_2O_2', 'propen2oxy+HO_2 \rightarrow acetone+O_2', 'HO_2+prod_7 \rightarrow well_5', 'propen1ol+HO_2 \rightarrow CH_2O+C_2H_3+H_2O_2', 'CH_3CHOH+HO_2 \rightarrow ethanol+O_2', 'ipropyl+HO_2 \rightarrow ipropyloxy+OH', 'ipropyloo+HO_2 \rightarrow ipropylooh+O_2', 'OH+HO_2 \rightarrow H_2O_2+O', 'HCO+HO_2 \rightarrow CH_2O+O_2', 'C_3H_6+HO_2 \rightarrow propen1yl+H_2O_2', 'allyl+HO_2 \rightarrow allyloxy+OH', 'allyl+HO_2 \rightarrow C_3H_6+O_2'};

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
ylim([10^-8, 10^0]);

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
figname = strcat(spe_name, '_sink_reaction_ratio_top_5', '.png');
print(fig, fullfile(file_dir, 'output', figname), '-r200', '-dpng');
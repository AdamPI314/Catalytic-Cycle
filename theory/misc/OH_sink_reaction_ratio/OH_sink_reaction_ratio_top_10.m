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
    0.3010    0.7450    0.9330 % 6th plot
    0.6350    0.0780    0.1840 % 7th plot
    0   0   1 % placeholder
    0   0.5   0 % placeholder
    0   0.75   0.75 % placeholder
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
OH_sink_R_idx = [1, 3, 4, 7, 21, 22, 29, 31, 32, 37, 39, 43, 44, 51, 52, 61, 62, 67, 79, 89, 90, 104, 111, 121, 122, 149, 152, 155, 161, 163, 164, 179, 183, 185, 186, 198, 209, 210, 212, 228, 236, 241, 243, 262, 267, 282, 290, 311, 312, 339, 361, 367, 369, 373, 377, 379, 383, 388, 411, 412, 424, 431, 447, 453, 461, 462, 464, 466, 471, 482, 501, 514, 517, 522, 524, 526, 533, 542, 544, 546, 567, 569, 571, 582, 589, 594, 599, 623, 628, 653, 654, 680, 685, 715, 717, 733, 735, 736, 738, 786, 815, 817, 819, 853, 862, 867, 884, 889, 898, 903, 907, 909, 1002, 1007, 1014, 1016, 1018, 1020, 1022, 1024, 1026, 1028, 1035, 1039, 1077, 1085, 1091, 1095, 1103, 1109, 1113, 1115, 1119, 1123, 1131, 1133, 1135, 1151, 1153, 1155, 1163, 1167, 1173, 1175, 1177, 1187, 1189, 1191, 1199, 1201, 1203, 1211, 1213, 1215, 1217, 1219, 1221];
for idx=1:length(OH_sink_R_idx)
    OH_sink_R_idx(idx) = OH_sink_R_idx(idx) + 1;
end

R_name = {'O+OH \rightarrow H+O_2', 'H+OH \rightarrow O+H_2', 'H_2+OH \rightarrow H_2O+H', 'OH+OH \rightarrow O+H_2O', 'OH+M \rightarrow O+H+M', 'H+OH+M \rightarrow H_2O+M', 'OH+OH \rightarrow HO_2+H', 'O_2+OH \rightarrow HO_2+O', 'HO_2+OH \rightarrow H_2O+O_2', 'OH+OH(+M) \rightarrow H_2O_2(+M)', 'H_2O+OH \rightarrow H_2O_2+H', 'OH+HO_2 \rightarrow H_2O_2+O', 'H_2O_2+OH \rightarrow HO_2+H_2O', 'CO_2+OH \rightarrow CO+HO_2', 'CO+OH \rightarrow CO_2+H', 'CO+OH \rightarrow HCO+O', 'HCO+OH \rightarrow CO+H_2O', 'CO_2+OH+H \rightarrow HCO+HO_2', 'formyloxy+OH \rightarrow formylooh', 'HCO+OH \rightarrow CH_2O+O', 'CH_2O+OH \rightarrow HCO+H_2O', 'HOCH_2O+OH \rightarrow HOCH_2OOH', 'CH_3O+OH \rightarrow CH_3+HO_2', 'CH_3+OH \rightarrow CH_4+O', 'CH_4+OH \rightarrow CH_3+H_2O', 'CH_3O+OH \rightarrow CH_3OO+H', 'CH_3OO+OH \rightarrow CH_3OH+O_2', 'CH_3O+OH \rightarrow CH_3OOH', 'CH_3+OH \rightarrow CH_2OH+H', 'CH_2O+OH \rightarrow CH_2OH+O', 'CH_2OH+OH \rightarrow CH_2O+H_2O', 'HOCH_2O+OH \rightarrow CH_2OH+HO_2', 'CH_3+OH \rightarrow CH_3O+H', 'CH_2O+OH \rightarrow CH_3O+O', 'CH_3O+OH \rightarrow CH_2O+H_2O', 'OH+CH_3(+M) \rightarrow CH_3OH(+M)', 'CH_2OH+OH \rightarrow CH_3OH+O', 'CH_3OH+OH \rightarrow CH_3O+H_2O', 'CH_3OH+OH \rightarrow CH_2OH+H_2O', 'CH_3+OH \rightarrow CH_2+H_2O', 'CH_2+OH \rightarrow CH_2O+H', 'HCO+OH \rightarrow CH_2+O_2', 'CH_2O+OH \rightarrow CH_2+HO_2', 'CH_2(S)+OH \rightarrow CH_2O+H', 'H+OH+CO \rightarrow CH_2(S)+O_2', 'CH_2+OH \rightarrow CH+H_2O', 'CH+OH \rightarrow HCO+H', 'C_2H_5+OH \rightarrow C_2H_6+O', 'C_2H_6+OH \rightarrow C_2H_5+H_2O', 'ethoxy+OH \rightarrow C_2H_5+HO_2', 'ethoxy+OH \rightarrow CH_3CH_2OOH', 'oxirane+OH \rightarrow C_2H_5+O_2', 'acetaldehyde+OH \rightarrow C_2H_5+O_2', 'acetaldehyde+OH \rightarrow CH_3CH_2OO', 'oxirane+OH \rightarrow CH_3CH_2OO', 'oxirane+OH \rightarrow CH_2CH_2OOH', 'acetaldehyde+OH \rightarrow CH_2CH_2OOH', 'oxirane+OH \rightarrow oxiranyl+H_2O', 'acetyl+OH \rightarrow acetaldehyde+O', 'acetaldehyde+OH \rightarrow acetyl+H_2O', 'acetaldehyde+OH \rightarrow vinoxy+H_2O', 'ketene+OH \rightarrow acetyl+O', 'acetyloxy+OH \rightarrow CH_3CO3H', 'CH_2O+CO+OH \rightarrow vinoxy+O_2', 'HCCO+OH \rightarrow ketene+O', 'ketene+OH \rightarrow HCCO+H_2O', 'ketene+OH \rightarrow CH_2OH+CO', 'HCCO+OH \rightarrow H_2+CO+CO', 'OH+CO+CO \rightarrow HCCO+O_2', 'C_2H_4+OH \rightarrow C_2H_3+H_2O', 'oxirane+OH \rightarrow C_2H_4+HO_2', 'C_2H_3+OH \rightarrow C_2H_2+H_2O', 'HCCO+OH \rightarrow C_2H_2+O_2', 'C_2H_2+OH \rightarrow ketene+H', 'C_2H_2+OH \rightarrow CH_3+CO', 'OH+C_2H_2 \rightarrow H+ethynol', 'C_2H_5+OH(+M) \rightarrow ethanol(+M)', 'ethanol+OH \rightarrow CH_2CH_2OH+H_2O', 'ethanol+OH \rightarrow CH_3CHOH+H_2O', 'ethanol+OH \rightarrow ethoxy+H_2O', 'CH_2CH_2OH+OH \rightarrow ethanol+O', 'CH_3CHOH+OH \rightarrow ethanol+O', 'ethoxy+OH \rightarrow ethanol+O', 'C_2H_4+OH \rightarrow CH_2CH_2OH', 'OH+CH_2O+CH_2O \rightarrow O_2C_2H_4OH', 'acetone+OH \rightarrow propen2oxy+H_2O', 'propen2oxy+OH \rightarrow acetone+O', 'CH_2CHCO+OH \rightarrow acrolein+O', 'acrolein+OH \rightarrow CH_2CHCO+H_2O', 'propionyl+OH \rightarrow propanal+O', 'propanal+OH \rightarrow propionyl+H_2O', 'CH_3OCH_3+OH \rightarrow CH_3OCH_2+H_2O', 'CH_3OCH_2+OH \rightarrow CH_3OCH_3+O', 'npropyloxy+OH \rightarrow npropylooh', 'ipropyloxy+OH \rightarrow ipropylooh', 'ipropyl+OH \rightarrow C_3H_8+O', 'npropyl+OH \rightarrow C_3H_8+O', 'C_3H_8+OH \rightarrow npropyl+H_2O', 'C_3H_8+OH \rightarrow ipropyl+H_2O', 'ipropyl+OH \rightarrow C_3H_6+H_2O', 'allyl+OH \rightarrow C_3H_6+O', 'propen1yl+OH \rightarrow C_3H_6+O', 'propen2yl+OH \rightarrow C_3H_6+O', 'propen1ol+OH \rightarrow C_3H_6+HO_2', 'allyl+OH \rightarrow acrolein+H+H', 'acrolein+OH \rightarrow allyl+O_2', 'propen1yl+OH \rightarrow C_2H_4+HCO+H', 'C_2H_4+HCO+OH \rightarrow propen1yl+HO_2', 'propen2yl+OH \rightarrow CH_3+ketene+H', 'CH_3+ketene+OH \rightarrow propen2yl+HO_2', 'npropyloxy+OH \rightarrow npropyl+HO_2', 'ipropyloxy+OH \rightarrow ipropyl+HO_2', 'propen1ol+OH \rightarrow CH_2O+C_2H_3+H_2O', 'CH_2O+C_2H_3+OH \rightarrow propen1ol+O', 'C_3H_6+OH \rightarrow allyl+H_2O', 'C_3H_6+OH \rightarrow propen1yl+H_2O', 'C_3H_6+OH \rightarrow propen2yl+H_2O', 'C_3H_6+OH \rightarrow allyl-alcohol+H', 'C_3H_6+OH \rightarrow ethenol+CH_3', 'C_3H_6+OH \rightarrow propen1ol+H', 'C_3H_6+OH \rightarrow propen2ol+H', 'C_3H_6+OH \rightarrow acetaldehyde+CH_3', 'allyloxy+OH \rightarrow allyl+HO_2', 'allyloxy+OH \rightarrow prod_2', 'OH+propoxide \rightarrow O_2+npropyl', 'OH+propoxide \rightarrow npropyloo', 'OH+propoxide \rightarrow QOOH_2', 'OH+propoxide \rightarrow QOOH_1', 'OH+propoxide \rightarrow O_2+ipropyl', 'OH+propoxide \rightarrow ipropyloo', 'OH+propoxide \rightarrow QOOH_3', 'OH+propoxide \rightarrow HO_2+C_3H_6', 'OH+OH+frag_1 \rightarrow O_2+QOOH_1', 'OH+prod_3 \rightarrow O_2+QOOH_1', 'OH+prod_3 \rightarrow O_2+QOOH_2', 'OH+OH+frag_4 \rightarrow O_2+QOOH_2', 'OH+OH+frag_5 \rightarrow O_2+QOOH_2', 'OH+prod_3 \rightarrow O_2+QOOH_3', 'OH+OH+frag_4 \rightarrow O_2+QOOH_3', 'OH+OH+frag_5 \rightarrow O_2+QOOH_3', 'OH+prod_1 \rightarrow well_1', 'OH+prod_3 \rightarrow well_1', 'OH+prod_3 \rightarrow well_2', 'OH+prod_4 \rightarrow well_2', 'OH+prod_5 \rightarrow well_2', 'OH+prod_3 \rightarrow well_3', 'OH+prod_4 \rightarrow well_3', 'OH+prod_5 \rightarrow well_3', 'OH+prod_3 \rightarrow well_5', 'OH+prod_4 \rightarrow well_5', 'OH+prod_5 \rightarrow well_5', 'propen1oxy+OH \rightarrow prod_6', 'propen2oxy+OH \rightarrow prod_7', 'frag_1+OH \rightarrow prod_1', 'frag_3+OH \rightarrow prod_3', 'frag_4+OH \rightarrow prod_4', 'frag_5+OH \rightarrow prod_5'};

OH_sink_R_mat = reaction_R_mat(:, OH_sink_R_idx);

% sort by the reaction rates around 0.5 tau, idx == 3500
[B,I] = sort(OH_sink_R_mat(3500, :),'descend');

for idx=1:length(I)
    r_idx = I(idx);
    if idx == 1
        R_total = OH_sink_R_mat(:, r_idx);
    else
        R_total = (R_total + OH_sink_R_mat(:, r_idx));
    end
end

N = 10;

% graph handler
H = gobjects(N);

for idx=1:N
    r_idx = I(idx);
    H(idx) = semilogy(time_vec, (OH_sink_R_mat(:, r_idx)) ./ (R_total), ...
        'LineWidth', 2); hold on;
end



%% settings
set(gca,'GridLineStyle','--');
xlabel('Time (seconds)', 'FontSize', 20);
ylabel('Fraction', 'FontSize', 20);
ylim([10^-10, 10^0]);

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
leg_h = legend([H(1);H(2);H(3);H(4);H(5);H(6);H(7);H(8);H(9);H(10)], ...
    R_name{1,I(1)},R_name{1,I(2)},R_name{1,I(3)},R_name{1,I(4)},R_name{1,I(5)},...
    R_name{1,I(6)},R_name{1,I(7)},R_name{1,I(8)},R_name{1,I(9)},R_name{1,I(10)});
set(leg_h, 'FontSize', 12, 'Box', 'off');
set(leg_h, 'Location', 'South')


%% save to file
figname = strcat('OH_sink_reaction_ratio_top_10', '.png');
print(fig, fullfile(file_dir, 'output', figname), '-r200', '-dpng');
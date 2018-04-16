%% global settings
spe_name = 'OH';

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
R_idx = [0, 2, 5, 6, 20, 23, 28, 30, 33, 36, 38, 42, 45, 50, 53, 60, 63, 66, 78, 88, 91, 105, 110, 120, 123, 148, 153, 154, 160, 162, 165, 178, 182, 184, 187, 199, 208, 211, 213, 229, 237, 240, 242, 263, 266, 283, 291, 310, 313, 338, 360, 366, 368, 372, 376, 378, 382, 389, 410, 413, 425, 430, 446, 452, 460, 463, 465, 467, 470, 483, 500, 515, 516, 523, 525, 527, 532, 543, 545, 547, 566, 568, 570, 583, 588, 595, 598, 622, 629, 652, 655, 681, 684, 714, 716, 732, 734, 737, 739, 787, 814, 816, 818, 852, 863, 866, 885, 888, 899, 902, 906, 908, 1003, 1006, 1015, 1017, 1019, 1021, 1023, 1025, 1027, 1029, 1034, 1038, 1076, 1084, 1090, 1094, 1102, 1108, 1112, 1114, 1118, 1122, 1130, 1132, 1134, 1150, 1152, 1154, 1162, 1166, 1172, 1174, 1176, 1186, 1188, 1190, 1198, 1200, 1202, 1210, 1212, 1214, 1216, 1218, 1220];
for idx=1:length(R_idx)
    R_idx(idx) = R_idx(idx) + 1;
end

R_name = {'H+O_2 \rightarrow O+OH', 'O+H_2 \rightarrow H+OH', 'H_2O+H \rightarrow H_2+OH', 'O+H_2O \rightarrow OH+OH', 'O+H+M \rightarrow OH+M', 'H_2O+M \rightarrow H+OH+M', 'HO_2+H \rightarrow OH+OH', 'HO_2+O \rightarrow O_2+OH', 'H_2O+O_2 \rightarrow HO_2+OH', 'H_2O_2(+M) \rightarrow OH+OH(+M)', 'H_2O_2+H \rightarrow H_2O+OH', 'H_2O_2+O \rightarrow OH+HO_2', 'HO_2+H_2O \rightarrow H_2O_2+OH', 'CO+HO_2 \rightarrow CO_2+OH', 'CO_2+H \rightarrow CO+OH', 'HCO+O \rightarrow CO+OH', 'CO+H_2O \rightarrow HCO+OH', 'HCO+HO_2 \rightarrow CO_2+OH+H', 'formylooh \rightarrow formyloxy+OH', 'CH_2O+O \rightarrow HCO+OH', 'HCO+H_2O \rightarrow CH_2O+OH', 'HOCH_2OOH \rightarrow HOCH_2O+OH', 'CH_3+HO_2 \rightarrow CH_3O+OH', 'CH_4+O \rightarrow CH_3+OH', 'CH_3+H_2O \rightarrow CH_4+OH', 'CH_3OO+H \rightarrow CH_3O+OH', 'CH_3OH+O_2 \rightarrow CH_3OO+OH', 'CH_3OOH \rightarrow CH_3O+OH', 'CH_2OH+H \rightarrow CH_3+OH', 'CH_2OH+O \rightarrow CH_2O+OH', 'CH_2O+H_2O \rightarrow CH_2OH+OH', 'CH_2OH+HO_2 \rightarrow HOCH_2O+OH', 'CH_3O+H \rightarrow CH_3+OH', 'CH_3O+O \rightarrow CH_2O+OH', 'CH_2O+H_2O \rightarrow CH_3O+OH', 'CH_3OH(+M) \rightarrow OH+CH_3(+M)', 'CH_3OH+O \rightarrow CH_2OH+OH', 'CH_3O+H_2O \rightarrow CH_3OH+OH', 'CH_2OH+H_2O \rightarrow CH_3OH+OH', 'CH_2+H_2O \rightarrow CH_3+OH', 'CH_2O+H \rightarrow CH_2+OH', 'CH_2+O_2 \rightarrow HCO+OH', 'CH_2+HO_2 \rightarrow CH_2O+OH', 'CH_2O+H \rightarrow CH_2(S)+OH', 'CH_2(S)+O_2 \rightarrow H+OH+CO', 'CH+H_2O \rightarrow CH_2+OH', 'HCO+H \rightarrow CH+OH', 'C_2H_6+O \rightarrow C_2H_5+OH', 'C_2H_5+H_2O \rightarrow C_2H_6+OH', 'C_2H_5+HO_2 \rightarrow ethoxy+OH', 'CH_3CH_2OOH \rightarrow ethoxy+OH', 'C_2H_5+O_2 \rightarrow oxirane+OH', 'C_2H_5+O_2 \rightarrow acetaldehyde+OH', 'CH_3CH_2OO \rightarrow acetaldehyde+OH', 'CH_3CH_2OO \rightarrow oxirane+OH', 'CH_2CH_2OOH \rightarrow oxirane+OH', 'CH_2CH_2OOH \rightarrow acetaldehyde+OH', 'oxiranyl+H_2O \rightarrow oxirane+OH', 'acetaldehyde+O \rightarrow acetyl+OH', 'acetyl+H_2O \rightarrow acetaldehyde+OH', 'vinoxy+H_2O \rightarrow acetaldehyde+OH', 'acetyl+O \rightarrow ketene+OH', 'CH_3CO_3H \rightarrow acetyloxy+OH', 'vinoxy+O_2 \rightarrow CH_2O+CO+OH', 'ketene+O \rightarrow HCCO+OH', 'HCCO+H_2O \rightarrow ketene+OH', 'CH_2OH+CO \rightarrow ketene+OH', 'H_2+CO+CO \rightarrow HCCO+OH', 'HCCO+O_2 \rightarrow OH+CO+CO', 'C_2H_3+H_2O \rightarrow C_2H_4+OH', 'C_2H_4+HO_2 \rightarrow oxirane+OH', 'C_2H_2+H_2O \rightarrow C_2H_3+OH', 'C_2H_2+O_2 \rightarrow HCCO+OH', 'ketene+H \rightarrow C_2H_2+OH', 'CH_3+CO \rightarrow C_2H_2+OH', 'H+ethynol \rightarrow OH+C_2H_2', 'ethanol(+M) \rightarrow C_2H_5+OH(+M)', 'CH_2CH_2OH+H_2O \rightarrow ethanol+OH', 'CH_3CHOH+H_2O \rightarrow ethanol+OH', 'ethoxy+H_2O \rightarrow ethanol+OH', 'ethanol+O \rightarrow CH_2CH_2OH+OH', 'ethanol+O \rightarrow CH_3CHOH+OH', 'ethanol+O \rightarrow ethoxy+OH', 'CH_2CH_2OH \rightarrow C_2H_4+OH', 'O_2C_2H_4OH \rightarrow OH+CH_2O+CH_2O', 'propen2oxy+H_2O \rightarrow acetone+OH', 'acetone+O \rightarrow propen2oxy+OH', 'acrolein+O \rightarrow CH_2CHCO+OH', 'CH_2CHCO+H_2O \rightarrow acrolein+OH', 'propanal+O \rightarrow propionyl+OH', 'propionyl+H_2O \rightarrow propanal+OH', 'CH_3OCH_2+H_2O \rightarrow CH_3OCH_3+OH', 'CH_3OCH_3+O \rightarrow CH_3OCH_2+OH', 'npropylooh \rightarrow npropyloxy+OH', 'ipropylooh \rightarrow ipropyloxy+OH', 'C_3H_8+O \rightarrow ipropyl+OH', 'C_3H_8+O \rightarrow npropyl+OH', 'npropyl+H_2O \rightarrow C_3H_8+OH', 'ipropyl+H_2O \rightarrow C_3H_8+OH', 'C_3H_6+H_2O \rightarrow ipropyl+OH', 'C_3H_6+O \rightarrow allyl+OH', 'C_3H_6+O \rightarrow propen1yl+OH', 'C_3H_6+O \rightarrow propen2yl+OH', 'C_3H_6+HO_2 \rightarrow propen1ol+OH', 'acrolein+H+H \rightarrow allyl+OH', 'allyl+O_2 \rightarrow acrolein+OH', 'C_2H_4+HCO+H \rightarrow propen1yl+OH', 'propen1yl+HO_2 \rightarrow C_2H_4+HCO+OH', 'CH_3+ketene+H \rightarrow propen2yl+OH', 'propen2yl+HO_2 \rightarrow CH_3+ketene+OH', 'npropyl+HO_2 \rightarrow npropyloxy+OH', 'ipropyl+HO_2 \rightarrow ipropyloxy+OH', 'CH_2O+C_2H_3+H_2O \rightarrow propen1ol+OH', 'propen1ol+O \rightarrow CH_2O+C_2H_3+OH', 'allyl+H_2O \rightarrow C_3H_6+OH', 'propen1yl+H_2O \rightarrow C_3H_6+OH', 'propen2yl+H_2O \rightarrow C_3H_6+OH', 'allyl-alcohol+H \rightarrow C_3H_6+OH', 'ethenol+CH_3 \rightarrow C_3H_6+OH', 'propen1ol+H \rightarrow C_3H_6+OH', 'propen2ol+H \rightarrow C_3H_6+OH', 'acetaldehyde+CH_3 \rightarrow C_3H_6+OH', 'allyl+HO_2 \rightarrow allyloxy+OH', 'prod_2 \rightarrow allyloxy+OH', 'O_2+npropyl \rightarrow OH+propoxide', 'npropyloo \rightarrow OH+propoxide', 'QOOH_2 \rightarrow OH+propoxide', 'QOOH_1 \rightarrow OH+propoxide', 'O_2+ipropyl \rightarrow OH+propoxide', 'ipropyloo \rightarrow OH+propoxide', 'QOOH_3 \rightarrow OH+propoxide', 'HO_2+C_3H_6 \rightarrow OH+propoxide', 'O_2+QOOH_1 \rightarrow OH+OH+frag_1', 'O_2+QOOH_1 \rightarrow OH+prod_3', 'O_2+QOOH_2 \rightarrow OH+prod_3', 'O_2+QOOH_2 \rightarrow OH+OH+frag_4', 'O_2+QOOH_2 \rightarrow OH+OH+frag_5', 'O_2+QOOH_3 \rightarrow OH+prod_3', 'O_2+QOOH_3 \rightarrow OH+OH+frag_4', 'O_2+QOOH_3 \rightarrow OH+OH+frag_5', 'well_1 \rightarrow OH+prod_1', 'well_1 \rightarrow OH+prod_3', 'well_2 \rightarrow OH+prod_3', 'well_2 \rightarrow OH+prod_4', 'well_2 \rightarrow OH+prod_5', 'well_3 \rightarrow OH+prod_3', 'well_3 \rightarrow OH+prod_4', 'well_3 \rightarrow OH+prod_5', 'well_5 \rightarrow OH+prod_3', 'well_5 \rightarrow OH+prod_4', 'well_5 \rightarrow OH+prod_5', 'prod_6 \rightarrow propen1oxy+OH', 'prod_7 \rightarrow propen2oxy+OH', 'prod_1 \rightarrow frag_1+OH', 'prod_3 \rightarrow frag_3+OH', 'prod_4 \rightarrow frag_4+OH', 'prod_5 \rightarrow frag_5+OH'};

R_mat = reaction_R_mat(:, R_idx);

% sort by the reaction rates around 0.5 tau, idx == 3550 for example
sort_axis = round(0.45 * length(time_vec));
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
    y_data = (R_mat(:, r_idx)) ./ (R_total);
    H(idx) = semilogy(time_vec(10:1:end), y_data(10:1:end), ...
        'LineWidth', 2); hold on;
end

%% settings
set(gca,'GridLineStyle','--');
xlabel('Time (seconds)', 'FontSize', 20);
ylabel('Fraction', 'FontSize', 20);
ylim([10^-2, 10^-0.1]);

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
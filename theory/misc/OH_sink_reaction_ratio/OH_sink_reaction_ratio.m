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
    1   0   0 % 8th plot, blue
    0   0   1 % placeholder
    0   0.5   0 % placeholder
    0   0.75   0.75 % placeholder
    0.7500         0    0.7500 % placeholder
    0.7500    0.7500         0 % placeholder
    0.2500    0.2500    0.2500 % placeholder
    0   0   1 % blue
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

OH_sink_R_name = {'O+OH=>H+O2', 'H+OH=>O+H2', 'H2+OH=>H2O+H', 'OH+OH=>O+H2O', 'OH+M=>O+H+M', 'H+OH+M=>H2O+M', 'OH+OH=>HO2+H', 'O2+OH=>HO2+O', 'HO2+OH=>H2O+O2', 'OH+OH(+M)=>H2O2(+M)', 'H2O+OH=>H2O2+H', 'OH+HO2=>H2O2+O', 'H2O2+OH=>HO2+H2O', 'CO2+OH=>CO+HO2', 'CO+OH=>CO2+H', 'CO+OH=>HCO+O', 'HCO+OH=>CO+H2O', 'CO2+OH+H=>HCO+HO2', 'formyloxy+OH=>formylooh', 'HCO+OH=>CH2O+O', 'CH2O+OH=>HCO+H2O', 'HOCH2O+OH=>HOCH2OOH', 'CH3O+OH=>CH3+HO2', 'CH3+OH=>CH4+O', 'CH4+OH=>CH3+H2O', 'CH3O+OH=>CH3OO+H', 'CH3OO+OH=>CH3OH+O2', 'CH3O+OH=>CH3OOH', 'CH3+OH=>CH2OH+H', 'CH2O+OH=>CH2OH+O', 'CH2OH+OH=>CH2O+H2O', 'HOCH2O+OH=>CH2OH+HO2', 'CH3+OH=>CH3O+H', 'CH2O+OH=>CH3O+O', 'CH3O+OH=>CH2O+H2O', 'OH+CH3(+M)=>CH3OH(+M)', 'CH2OH+OH=>CH3OH+O', 'CH3OH+OH=>CH3O+H2O', 'CH3OH+OH=>CH2OH+H2O', 'CH3+OH=>CH2+H2O', 'CH2+OH=>CH2O+H', 'HCO+OH=>CH2+O2', 'CH2O+OH=>CH2+HO2', 'CH2(S)+OH=>CH2O+H', 'H+OH+CO=>CH2(S)+O2', 'CH2+OH=>CH+H2O', 'CH+OH=>HCO+H', 'C2H5+OH=>C2H6+O', 'C2H6+OH=>C2H5+H2O', 'ethoxy+OH=>C2H5+HO2', 'ethoxy+OH=>CH3CH2OOH', 'oxirane+OH=>C2H5+O2', 'acetaldehyde+OH=>C2H5+O2', 'acetaldehyde+OH=>CH3CH2OO', 'oxirane+OH=>CH3CH2OO', 'oxirane+OH=>CH2CH2OOH', 'acetaldehyde+OH=>CH2CH2OOH', 'oxirane+OH=>oxiranyl+H2O', 'acetyl+OH=>acetaldehyde+O', 'acetaldehyde+OH=>acetyl+H2O', 'acetaldehyde+OH=>vinoxy+H2O', 'ketene+OH=>acetyl+O', 'acetyloxy+OH=>CH3CO3H', 'CH2O+CO+OH=>vinoxy+O2', 'HCCO+OH=>ketene+O', 'ketene+OH=>HCCO+H2O', 'ketene+OH=>CH2OH+CO', 'HCCO+OH=>H2+CO+CO', 'OH+CO+CO=>HCCO+O2', 'C2H4+OH=>C2H3+H2O', 'oxirane+OH=>C2H4+HO2', 'C2H3+OH=>C2H2+H2O', 'HCCO+OH=>C2H2+O2', 'C2H2+OH=>ketene+H', 'C2H2+OH=>CH3+CO', 'OH+C2H2=>H+ethynol', 'C2H5+OH(+M)=>ethanol(+M)', 'ethanol+OH=>CH2CH2OH+H2O', 'ethanol+OH=>CH3CHOH+H2O', 'ethanol+OH=>ethoxy+H2O', 'CH2CH2OH+OH=>ethanol+O', 'CH3CHOH+OH=>ethanol+O', 'ethoxy+OH=>ethanol+O', 'C2H4+OH=>CH2CH2OH', 'OH+CH2O+CH2O=>O2C2H4OH', 'acetone+OH=>propen2oxy+H2O', 'propen2oxy+OH=>acetone+O', 'CH2CHCO+OH=>acrolein+O', 'acrolein+OH=>CH2CHCO+H2O', 'propionyl+OH=>propanal+O', 'propanal+OH=>propionyl+H2O', 'CH3OCH3+OH=>CH3OCH2+H2O', 'CH3OCH2+OH=>CH3OCH3+O', 'npropyloxy+OH=>npropylooh', 'ipropyloxy+OH=>ipropylooh', 'ipropyl+OH=>C3H8+O', 'npropyl+OH=>C3H8+O', 'C3H8+OH=>npropyl+H2O', 'C3H8+OH=>ipropyl+H2O', 'ipropyl+OH=>C3H6+H2O', 'allyl+OH=>C3H6+O', 'propen1yl+OH=>C3H6+O', 'propen2yl+OH=>C3H6+O', 'propen1ol+OH=>C3H6+HO2', 'allyl+OH=>acrolein+H+H', 'acrolein+OH=>allyl+O2', 'propen1yl+OH=>C2H4+HCO+H', 'C2H4+HCO+OH=>propen1yl+HO2', 'propen2yl+OH=>CH3+ketene+H', 'CH3+ketene+OH=>propen2yl+HO2', 'npropyloxy+OH=>npropyl+HO2', 'ipropyloxy+OH=>ipropyl+HO2', 'propen1ol+OH=>CH2O+C2H3+H2O', 'CH2O+C2H3+OH=>propen1ol+O', 'C3H6+OH=>allyl+H2O', 'C3H6+OH=>propen1yl+H2O', 'C3H6+OH=>propen2yl+H2O', 'C3H6+OH=>allyl-alcohol+H', 'C3H6+OH=>ethenol+CH3', 'C3H6+OH=>propen1ol+H', 'C3H6+OH=>propen2ol+H', 'C3H6+OH=>acetaldehyde+CH3', 'allyloxy+OH=>allyl+HO2', 'allyloxy+OH=>prod_2', 'OH+propoxide=>O2+npropyl', 'OH+propoxide=>npropyloo', 'OH+propoxide=>QOOH_2', 'OH+propoxide=>QOOH_1', 'OH+propoxide=>O2+ipropyl', 'OH+propoxide=>ipropyloo', 'OH+propoxide=>QOOH_3', 'OH+propoxide=>HO2+C3H6', 'OH+OH+frag_1=>O2+QOOH_1', 'OH+prod_3=>O2+QOOH_1', 'OH+prod_3=>O2+QOOH_2', 'OH+OH+frag_4=>O2+QOOH_2', 'OH+OH+frag_5=>O2+QOOH_2', 'OH+prod_3=>O2+QOOH_3', 'OH+OH+frag_4=>O2+QOOH_3', 'OH+OH+frag_5=>O2+QOOH_3', 'OH+prod_1=>well_1', 'OH+prod_3=>well_1', 'OH+prod_3=>well_2', 'OH+prod_4=>well_2', 'OH+prod_5=>well_2', 'OH+prod_3=>well_3', 'OH+prod_4=>well_3', 'OH+prod_5=>well_3', 'OH+prod_3=>well_5', 'OH+prod_4=>well_5', 'OH+prod_5=>well_5', 'propen1oxy+OH=>prod_6', 'propen2oxy+OH=>prod_7', 'frag_1+OH=>prod_1', 'frag_3+OH=>prod_3', 'frag_4+OH=>prod_4', 'frag_5+OH=>prod_5'};

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

% graph handler
H = gobjects(length(I));

for idx=1:length(I)
% for idx=1:10
    r_idx = I(idx);
    H(idx) = semilogy(time_vec, (OH_sink_R_mat(:, r_idx)) ./ (R_total), ...
        'LineWidth', 1.25); hold on;
end



%% settings
set(gca,'GridLineStyle','--');
xlabel('Time (seconds)', 'FontSize', 20);
ylabel('Fraction', 'FontSize', 20);
ylim([10^-60, 10^0]);

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
% leg_h = legend([p1; p2; p3; p4; p5; p6],'[nROO]/[nR]','K_{eq}','[nROO]/[QOOH_1]','K_{eq}','[O_2QOOH_1]/[QOOH_1]','K_{eq}');
% set(leg_h, 'FontSize', 14, 'Box', 'off');
% set(leg_h, 'Location', 'South')


%% save to file
figname = strcat('OH_sink_reaction_ratio_all', '.png');
print(fig, fullfile(file_dir, 'output', figname), '-r200', '-dpng');
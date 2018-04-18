%% global settings
spe_name = 'H2O2';

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
R_idx = [69, 97, 107, 109, 111, 113, 115, 117, 118, 120, 122, 124, 127, 129, 133, 136, 141, 160, 182, 192, 199, 221, 225, 226, 229, 231, 232, 238, 264, 272, 317, 333, 345, 352, 384, 399, 407, 417, 426, 433, 440, 448, 454, 478, 485, 511, 524, 530, 573, 575, 577, 592, 601, 610, 635, 657, 678, 691, 698, 718, 745, 747, 784, 788, 790, 803, 810, 826, 841, 843, 845, 875, 881, 893, 896, 898, 902, 938, 940, 973, 983, 997, 1013, 1022, 1028, 1226];
for idx=1:length(R_idx)
    R_idx(idx) = R_idx(idx) + 1;
end

R_name = {'CO+CH_4=>HCO+CH_3', 'HCO+CH_4=>CH_2O+CH_3', 'CH_2O+H=>CH_3+O', 'CH_3O+O=>CH_3+O_2', 'CH_3O+OH=>CH_3+HO_2', 'CH_4+O_2=>CH_3+HO_2', 'C_2H_6(+M)=>CH_3+CH_3(+M)', 'CH_4(+M)=>CH_3+H(+M)', 'CH_4+H=>CH_3+H_2', 'CH_4+O=>CH_3+OH', 'CH_4+OH=>CH_3+H_2O', 'CH_4+HO_2=>CH_3+H_2O_2', 'CH_4+CH_3O=>CH_3+CH_3OH', 'CH_2O+CH_4=>CH_3O+CH_3', 'CH_3OO(+M)=>CH_3+O_2(+M)', 'CH_4+CH_3OO=>CH_3+CH_3OOH', 'CH_3O+CH_3O=>CH_3OO+CH_3', 'CH_2OH+H=>CH_3+OH', 'CH_3O+H=>CH_3+OH', 'CH_3O+CO=>CH_3+CO_2', 'CH_3OH(+M)=>OH+CH_3(+M)', 'CH_2OH+CH_4=>CH_3OH+CH_3', 'H+C_2H_5=>CH_3+CH_3', 'CH_4+CH_2=>CH_3+CH_3', 'CH_2+H_2O=>CH_3+OH', 'C_2H_4+H=>CH_3+CH_2', 'CH_2+H(+M)=>CH_3(+M)', 'CH_2+H_2=>H+CH_3', 'CH_2(S)+H_2=>CH_3+H', 'CH_2(S)+C_2H_6=>CH_3+C_2H_5', 'C_2H_5+CH_4=>C_2H_6+CH_3', 'CH_4+C_2H_4=>CH_3+C_2H_5', 'ethoxy=>CH_3+CH_2O', 'CH_4+CH_3CH_2OO=>CH_3+CH_3CH_2OOH', 'oxirane=>CH_3+HCO', 'oxiranyl+CH_4=>oxirane+CH_3', 'acetaldehyde=>CH_3+HCO', 'acetyl+CH_4=>acetaldehyde+CH_3', 'acetyl(+M)=>CH_3+CO(+M)', 'ketene+CH_4=>acetyl+CH_3', 'CH_4+acetylperoxy=>CH_3+CH_3CO_3H', 'acetyloxy+M=>CH_3+CO_2+M', 'ketene+H=>CH_3+CO', 'C_2H_4+O=>CH_3+HCO', 'C_2H_3+CH_4=>C_2H_4+CH_3', 'CH_4+C_2H_2=>CH_3+C_2H_3', 'C_2H_2+OH=>CH_3+CO', 'ethanol(+M)=>CH_2OH+CH_3(+M)', 'CH_2CH_2OH+CH_4=>ethanol+CH_3', 'CH_3CHOH+CH_4=>ethanol+CH_3', 'ethoxy+CH_4=>ethanol+CH_3', 'acetone(+M)=>acetyl+CH_3(+M)', 'propen2oxy+CH_4=>acetone+CH_3', 'propen2oxy=>ketene+CH_3', 'CH_2CHCO+CH_4=>acrolein+CH_3', 'propionyl+CH_4=>propanal+CH_3', 'CH_3OCH_3(+M)=>CH_3+CH_3O(+M)', 'CH_3OCH_2+CH_4=>CH_3OCH_3+CH_3', 'CH_3OCH_2=>CH_2O+CH_3', 'C_3H_8(+M)=>CH_3+C_2H_5(+M)', 'CH_4+ipropyl=>CH_3+C_3H_8', 'CH_4+npropyl=>CH_3+C_3H_8', 'ipropyl+H=>C_2H_5+CH_3', 'ipropyl+O=>acetaldehyde+CH_3', 'npropyl=>CH_3+C_2H_4', 'C_3H_6(+M)=>C_2H_3+CH_3(+M)', 'C_3H_6+O=>ketene+CH_3+H', 'C_3H_6+H=>C_2H_4+CH_3', 'allyl+CH_4=>C_3H_6+CH_3', 'propen1yl+CH_4=>C_3H_6+CH_3', 'propen2yl+CH_4=>C_3H_6+CH_3', 'allyl=>C_2H_2+CH_3', 'propen1yl=>C_2H_2+CH_3', 'propen2yl=>C_2H_2+CH_3', 'propen2yl+O=>CH_3+ketene', 'propen2yl+OH=>CH_3+ketene+H', 'propen2yl+HO_2=>CH_3+ketene+OH', 'CH_4+npropyloo=>CH_3+npropylooh', 'CH_4+ipropyloo=>CH_3+ipropylooh', 'ipropyloxy+CH_3O=>ipropyloo+CH_3', 'npropyloxy+CH_3O=>npropyloo+CH_3', 'ipropyloxy=>CH_3+acetaldehyde', 'CH_2O+C_2H_3+CH_4=>propen1ol+CH_3', 'C_3H_6+OH=>ethenol+CH_3', 'C_3H_6+OH=>acetaldehyde+CH_3', 'frag_5=>CH_3+glyoxal'};

%% chattering reactions
chattering_R_idx = [1068, 1069, 1096, 1097, 1080, 1081, 132, 133, 348, 349];
R_idx_reduced = {};
R_name_reduced = {};
for i=1:length(R_idx)
    if ~ ismember(R_idx(i)-1, chattering_R_idx)
        R_idx_reduced{end+1} = R_idx(i);
        R_name_reduced{end+1} = R_name{1,i};
    end
end
R_idx = cell2mat(R_idx_reduced);
R_name = R_name_reduced;

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
ylim([10^-5, 10^0]);

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
figname = strcat(spe_name, '_formation_pathway_prob_vs_time', '.png');
print(fig, fullfile(pic_dir, figname), '-r200', '-dpng');
%% global settings
spe_name = 'CH3';

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
R_idx = [68, 96, 106, 108, 110, 112, 114, 116, 119, 121, 123, 125, 126, 128, 132, 137, 140, 161, 183, 193, 198, 220, 224, 227, 228, 230, 233, 239, 265, 273, 316, 332, 344, 353, 385, 398, 406, 416, 427, 432, 441, 449, 455, 479, 484, 510, 525, 531, 572, 574, 576, 593, 600, 611, 634, 656, 679, 690, 699, 719, 744, 746, 785, 789, 791, 802, 811, 827, 840, 842, 844, 874, 880, 892, 897, 899, 903, 939, 941, 972, 982, 996, 1012, 1023, 1029, 1227];
for idx=1:length(R_idx)
    R_idx(idx) = R_idx(idx) + 1;
end

R_name = {'HCO+CH_3 \rightarrow CO+CH_4', 'CH_2O+CH_3 \rightarrow HCO+CH_4', 'CH_3+O \rightarrow CH_2O+H', 'CH_3+O_2 \rightarrow CH_3O+O', 'CH_3+HO_2 \rightarrow CH_3O+OH', 'CH_3+HO_2 \rightarrow CH_4+O_2', 'CH_3+CH_3(+M) \rightarrow C_2H_6(+M)', 'CH_3+H(+M) \rightarrow CH_4(+M)', 'CH_3+H_2 \rightarrow CH_4+H', 'CH_3+OH \rightarrow CH_4+O', 'CH_3+H_2O \rightarrow CH_4+OH', 'CH_3+H_2O_2 \rightarrow CH_4+HO_2', 'CH_3+CH_3OH \rightarrow CH_4+CH_3O', 'CH_3O+CH_3 \rightarrow CH_2O+CH_4', 'CH_3+O_2(+M) \rightarrow CH_3OO(+M)', 'CH_3+CH_3OOH \rightarrow CH_4+CH_3OO', 'CH_3OO+CH_3 \rightarrow CH_3O+CH_3O', 'CH_3+OH \rightarrow CH_2OH+H', 'CH_3+OH \rightarrow CH_3O+H', 'CH_3+CO_2 \rightarrow CH_3O+CO', 'OH+CH_3(+M) \rightarrow CH_3OH(+M)', 'CH_3OH+CH_3 \rightarrow CH_2OH+CH_4', 'CH_3+CH_3 \rightarrow H+C_2H_5', 'CH_3+CH_3 \rightarrow CH_4+CH_2', 'CH_3+OH \rightarrow CH_2+H_2O', 'CH_3+CH_2 \rightarrow C_2H_4+H', 'CH_3(+M) \rightarrow CH_2+H(+M)', 'H+CH_3 \rightarrow CH_2+H_2', 'CH_3+H \rightarrow CH_2(S)+H_2', 'CH_3+C_2H_5 \rightarrow CH_2(S)+C_2H_6', 'C_2H_6+CH_3 \rightarrow C_2H_5+CH_4', 'CH_3+C_2H_5 \rightarrow CH_4+C_2H_4', 'CH_3+CH_2O \rightarrow ethoxy', 'CH_3+CH_3CH_2OOH \rightarrow CH_4+CH_3CH_2OO', 'CH_3+HCO \rightarrow oxirane', 'oxirane+CH_3 \rightarrow oxiranyl+CH_4', 'CH_3+HCO \rightarrow acetaldehyde', 'acetaldehyde+CH_3 \rightarrow acetyl+CH_4', 'CH_3+CO(+M) \rightarrow acetyl(+M)', 'acetyl+CH_3 \rightarrow ketene+CH_4', 'CH_3+CH_3CO_3H \rightarrow CH_4+acetylperoxy', 'CH_3+CO_2+M \rightarrow acetyloxy+M', 'CH_3+CO \rightarrow ketene+H', 'CH_3+HCO \rightarrow C_2H_4+O', 'C_2H_4+CH_3 \rightarrow C_2H_3+CH_4', 'CH_3+C_2H_3 \rightarrow CH_4+C_2H_2', 'CH_3+CO \rightarrow C_2H_2+OH', 'CH_2OH+CH_3(+M) \rightarrow ethanol(+M)', 'ethanol+CH_3 \rightarrow CH_2CH_2OH+CH_4', 'ethanol+CH_3 \rightarrow CH_3CHOH+CH_4', 'ethanol+CH_3 \rightarrow ethoxy+CH_4', 'acetyl+CH_3(+M) \rightarrow acetone(+M)', 'acetone+CH_3 \rightarrow propen2oxy+CH_4', 'ketene+CH_3 \rightarrow propen2oxy', 'acrolein+CH_3 \rightarrow CH_2CHCO+CH_4', 'propanal+CH_3 \rightarrow propionyl+CH_4', 'CH_3+CH_3O(+M) \rightarrow CH_3OCH_3(+M)', 'CH_3OCH_3+CH_3 \rightarrow CH_3OCH_2+CH_4', 'CH_2O+CH_3 \rightarrow CH_3OCH_2', 'CH_3+C_2H_5(+M) \rightarrow C_3H_8(+M)', 'CH_3+C_3H_8 \rightarrow CH_4+ipropyl', 'CH_3+C_3H_8 \rightarrow CH_4+npropyl', 'C_2H_5+CH_3 \rightarrow ipropyl+H', 'acetaldehyde+CH_3 \rightarrow ipropyl+O', 'CH_3+C_2H_4 \rightarrow npropyl', 'C_2H_3+CH_3(+M) \rightarrow C_3H_6(+M)', 'ketene+CH_3+H \rightarrow C_3H_6+O', 'C_2H_4+CH_3 \rightarrow C_3H_6+H', 'C_3H_6+CH_3 \rightarrow allyl+CH_4', 'C_3H_6+CH_3 \rightarrow propen1yl+CH_4', 'C_3H_6+CH_3 \rightarrow propen2yl+CH_4', 'C_2H_2+CH_3 \rightarrow allyl', 'C_2H_2+CH_3 \rightarrow propen1yl', 'C_2H_2+CH_3 \rightarrow propen2yl', 'CH_3+ketene \rightarrow propen2yl+O', 'CH_3+ketene+H \rightarrow propen2yl+OH', 'CH_3+ketene+OH \rightarrow propen2yl+HO_2', 'CH_3+npropylooh \rightarrow CH_4+npropyloo', 'CH_3+ipropylooh \rightarrow CH_4+ipropyloo', 'ipropyloo+CH_3 \rightarrow ipropyloxy+CH_3O', 'npropyloo+CH_3 \rightarrow npropyloxy+CH_3O', 'CH_3+acetaldehyde \rightarrow ipropyloxy', 'propen1ol+CH_3 \rightarrow CH_2O+C_2H_3+CH_4', 'ethenol+CH_3 \rightarrow C_3H_6+OH', 'acetaldehyde+CH_3 \rightarrow C_3H_6+OH', 'CH_3+glyoxal \rightarrow frag_5'};
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
sort_axis = round(0.1 * length(time_vec));
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
ylim([10^-3.5, 10^0]);

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
print(fig, fullfile(pic_dir, figname), '-r200', '-dpng');
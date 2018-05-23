%% global settings
spe_name = 'CH2O';

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
R_idx = [73, 76, 82, 84, 86, 88, 90, 92, 94, 96, 98, 107, 129, 131, 134, 145, 157, 159, 163, 165, 167, 169, 173, 175, 177, 181, 185, 187, 189, 191, 197, 217, 237, 243, 263, 271, 293, 298, 344, 350, 442, 453, 505, 589, 699, 701, 702, 865, 901, 914, 918, 992, 1001, 1003, 1005, 1007, 1009, 1011, 1013, 1041, 1051, 1059, 1064, 1066, 1223, 1225];
for idx=1:length(R_idx)
    R_idx(idx) = R_idx(idx) + 1;
end

R_name = {'CH_2O+CO \rightarrow HCO+HCO', 'CH_2O+formylperoxy \rightarrow HCO+formylooh', 'CH_2O+M \rightarrow HCO+H+M', 'CH_2O+M \rightarrow CO+H_2+M', 'CH_2O+H \rightarrow HCO+H_2', 'CH_2O+O \rightarrow HCO+OH', 'CH_2O+OH \rightarrow HCO+H_2O', 'CH_2O+O_2 \rightarrow HCO+HO_2', 'CH_2O+HO_2 \rightarrow HCO+H_2O_2', 'CH_2O+CH_3 \rightarrow HCO+CH_4', 'CH_2O+HO_2 \rightarrow OCH_2OOH', 'CH_2O+H \rightarrow CH_3+O', 'CH_2O+CH_4 \rightarrow CH_3O+CH_3', 'CH_2O+H_2 \rightarrow CH_3O+H', 'CH_3OO+CH_2O \rightarrow CH_3OOH+HCO', 'CH_2O+CH_3OH+O_2 \rightarrow CH_3OO+CH_3OO', 'CH_2O+H+M \rightarrow CH_2OH+M', 'CH_2O+H_2 \rightarrow CH_2OH+H', 'CH_2O+OH \rightarrow CH_2OH+O', 'CH_2O+H_2O \rightarrow CH_2OH+OH', 'CH_2O+HO_2 \rightarrow CH_2OH+O_2', 'CH_2O+H_2O_2 \rightarrow CH_2OH+HO_2', 'CH_2O+CH_2O \rightarrow CH_2OH+HCO', 'CH_3OH+CH_2O \rightarrow 2CH_2OH', 'CH_3OH+CH_2O \rightarrow CH_2OH+CH_3O', 'CH_2O+H+M \rightarrow CH_3O+M', 'CH_2O+OH \rightarrow CH_3O+O', 'CH_2O+H_2O \rightarrow CH_3O+OH', 'CH_2O+HO_2 \rightarrow CH_3O+O_2', 'CH_2O+H_2O_2 \rightarrow CH_3O+HO_2', 'CH_3OH+CH_2O \rightarrow 2CH_3O', 'CH_2OH+CH_2O \rightarrow CH_3OH+HCO', 'CH_2O+H \rightarrow CH_2+OH', 'CH_2O+OH \rightarrow CH_2+HO_2', 'CH_2O+H \rightarrow CH_2(S)+OH', 'CH_2O+CO \rightarrow CH_2(S)+CO_2', 'H+CH_2O \rightarrow CH+H_2O', 'CH+CH_2O \rightarrow H+ketene', 'CH_3+CH_2O \rightarrow ethoxy', 'CH_3CH_2OO+CH_2O \rightarrow CH_3CH_2OOH+HCO', 'CH_2O+acetylperoxy \rightarrow HCO+CH_3CO_3H', 'CH_2O+CO+OH \rightarrow vinoxy+O_2', 'HCO+CH_2O \rightarrow C_2H_3+O_2', 'OH+CH_2O+CH_2O \rightarrow O_2C_2H_4OH', 'CH_2O+CH_3 \rightarrow CH_3OCH_2', 'CH_3OCH_3+CH_2O \rightarrow CH_3OCH_2+CH_3O', 'CH_3OCH_2+CH_2O \rightarrow CH_3OCH_3+HCO', 'acetyl+CH_2O \rightarrow allyl+O_2', 'acetyl+CH_2O \rightarrow propen2yl+O_2', 'npropyloo+CH_2O \rightarrow npropylooh+HCO', 'ipropyloo+CH_2O \rightarrow ipropylooh+HCO', 'C_2H_5+CH_2O \rightarrow npropyloxy', 'C_2H_4+CH_2O \rightarrow propen1ol', 'CH_2O+C_2H_3+H_2O \rightarrow propen1ol+OH', 'CH_2O+C_2H_3+H_2 \rightarrow propen1ol+H', 'CH_2O+C_2H_3+OH \rightarrow propen1ol+O', 'CH_2O+C_2H_3+H_2O_2 \rightarrow propen1ol+HO_2', 'CH_2O+C_2H_3+CH_3OOH \rightarrow propen1ol+CH_3OO', 'CH_2O+C_2H_3+CH_4 \rightarrow propen1ol+CH_3', 'C_2H_3+CH_2O \rightarrow allyloxy', 'C_2H_3+CH_2O \rightarrow vinoxylmethyl', 'C_2H_3+CH_2O \rightarrow formylethyl', 'C_2H_3+CH_2O \rightarrow acrolein+H', 'C_2H_3+CH_2O \rightarrow C_2H_4+HCO', 'vinoxy+CH_2O \rightarrow frag_1', 'acetyl+CH_2O \rightarrow frag_4'};

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
ylim([10^-3, 10^0]);

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
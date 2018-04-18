%% global settings
spe_name = 'H2O2';
N_time = 25;
end_t = '0.9';
tau = 0.777660157519;
sort_time = '0.01';

%% Current file directory
file_dir = fullfile(fileparts(mfilename('fullpath')));
pic_dir = fullfile(fileparts(mfilename('fullpath')));

% import time
fn_time = fullfile(file_dir, 'pathway_time_candidate_S62_H_0.9_100000000_10000_100.csv');
delimiter = ',';
formatStr = '';
for i=1:N_time
    formatStr = strcat(formatStr, '%f');
end
formatStr = strcat(formatStr, '%[^\n\r]');
formatSpec = char(formatStr);
fileID = fopen(fn_time,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'ReturnOnError', false);
fclose(fileID);
time_vec = [dataArray{1:end-1}];
time_vec = time_vec(1, :);
%% Clear temporary variables
clearvars fn_time delimiter formatSpec fileID dataArray ans;

% import pathP
fn_pathP = fullfile(file_dir, 'pathway_prob_S62_H_0.9_100000000_10000_100.csv');
delimiter = ',';
formatStr = '';
for i=1:N_time
    formatStr = strcat(formatStr, '%f');
end
formatStr = strcat(formatStr, '%[^\n\r]');
formatSpec = char(formatStr);
fileID = fopen(fn_pathP,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'ReturnOnError', false);
fclose(fileID);
pathP_mat = [dataArray{1:end-1}];
%% Clear temporary variables
clearvars fn_pathP delimiter formatSpec fileID dataArray ans;

%% plot
% sort by the reaction rates around 0.5 tau, idx == 3550 for example
sort_axis = round(str2double(sort_time) * length(time_vec));
if sort_axis <= 0
    sort_axis = 1;
end
    
[B,I] = sort(pathP_mat(sort_axis, :),'descend');

N = 5;
% graph handler
H = gobjects(N);

for idx=1:N
    r_idx = I(idx);
    y_data = pathP_mat(r_idx, :);
    H(idx) = semilogy(time_vec * tau, y_data, ...
        'LineWidth', 2); hold on;
end

%% settings
set(gca,'GridLineStyle','--');
xlabel('Time (seconds)', 'FontSize', 20);
ylabel('Pathway Probability', 'FontSize', 20);
grid on;

% %% save to file
% figname = strcat(spe_name, '_formation_pathway_prob_vs_time', '.png');
% print(fig, fullfile(pic_dir, figname), '-r200', '-dpng');
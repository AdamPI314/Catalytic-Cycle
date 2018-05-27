%% global settings
file_dir = fullfile(fileparts(mfilename('fullpath')));
pic_dir = fullfile(fileparts(mfilename('fullpath')));

% marker
% markers = {'+' , 'o' , '*' , 'x' , 'square' , 'diamond' , 'v' , '^' , '>' , '<' , 'pentagram' , 'hexagram' , '.', 'none'};
markers = {'+' , 'o' , '*' , 'x' , 'square' , 'diamond' , 'v' , '^' , '>' , '<' , 'pentagram' , 'hexagram' , '.'};
% markers = {'None'};
prefix = '';
spe_idx = 14;
spe_name = 'CO';
tau = 0.777660157519;
% stage 1A end time
end_t = '0.25718313951098054';
% stage 1B end time
% end_t = '0.720108899222239';
n_path = 100;

%% import time
f_n_scc_time = fullfile(file_dir, [prefix, 'scc_time.csv']);
delimiter = {''};
formatSpec = '%f%[^\n\r]';
fileID = fopen(f_n_scc_time,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string',  'ReturnOnError', false);
fclose(fileID);
scc_time_vec = dataArray{:, 1};
clearvars f_n_scc_time delimiter formatSpec fileID dataArray ans;

%% import spe conc
f_n_conc = fullfile(file_dir, [prefix, 'spe_conc_converted_to_pp.csv']);
delimiter = ',';
formatStr = '';
for i=1:n_path
    formatStr = strcat(formatStr, '%f');
end
formatStr = strcat(formatStr, '%[^\n\r]');
formatSpec = char(formatStr);
fileID = fopen(f_n_conc,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN,  'ReturnOnError', false);
fclose(fileID);
spe_conc_in_pp = [dataArray{1:end-1}];
clearvars f_n_conc delimiter formatSpec fileID dataArray ans;

%% import pathway prob
f_n_pp = fullfile(file_dir, [prefix, 'scc_path_prob.csv']);
delimiter = ',';
formatStr = '';
for i=1:n_path
    formatStr = strcat(formatStr, '%f');
end
formatStr = strcat(formatStr, '%[^\n\r]');
formatSpec = char(formatStr);
fileID = fopen(f_n_pp,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN,  'ReturnOnError', false);
fclose(fileID);
scc_pp = [dataArray{1:end-1}];
clearvars f_n_pp delimiter formatSpec fileID dataArray ans;

%% global properties
idx_vec = [1 5 10 50];
N = length(idx_vec);
y_label_str = 'Normalized [X]';

%% plot
fig = figure();

%% plot
colors = lines(N + 1);
% colors = colorcube(N);

str_name = cell(N + 1,1);
for idx=1:N
    str_name{idx, 1} = ['N = ', num2str(idx_vec(idx))];
end
str_name{N+1, 1} = 'EXACT';

data_x = scc_time_vec * tau;
color_idx = 1;
for idx =1:N+1
    if idx == N+1
        data_y = spe_conc_in_pp(:, 14 + 1);
    else
        A = scc_pp(:, 1:idx_vec(idx));
        data_y = sum(A, 2);
    end

    plot(data_x, data_y, ...
    'LineWidth', 1, ...
    'color', colors(mod(color_idx-1, length(colors))+ 1, :), ...
    'marker', markers{1, mod(color_idx-1, length(markers))+ 1}, ...
    'HandleVisibility', 'on');
    hold on;
    color_idx = color_idx + 1;
end

%% settings
set(gca,'GridLineStyle','--');
xlabel('$t$ (seconds)', 'Interpreter','latex', 'FontSize', 20);
ylabel(y_label_str, 'Interpreter','latex', 'FontSize', 20);
% yticks([0 0.01 0.02 0.03 0.04 0.05]);
xlim([0,  str2double(end_t) * tau]);
% xlim([0,  str2double('0.25718313951098054') * tau]);
grid on;

%  legend
leg_h = legend(str_name, 'Interpreter','latex');
set(leg_h, 'FontSize', 12, 'Box', 'off');
l_pos = get(leg_h, 'Position');
set(leg_h, 'Position', [0 l_pos(2) 0.5 l_pos(4)]);
% set(leg_h, 'Location', 'West');

%% text
a_x = gca;
t_x = a_x.XLim(1) + 0.45*(a_x.XLim(2) - a_x.XLim(1));
t_y = a_x.YLim(1) + 0.80*(a_x.YLim(2) - a_x.YLim(1));
text(t_x, t_y, ['CO'], 'Interpreter','latex', 'FontSize', 20);

%% save to file
figname = strcat('spe_conc_path_prob_convergence_', end_t, '_', spe_name, '.png');
print(fig, fullfile(pic_dir, figname), '-r200', '-dpng');

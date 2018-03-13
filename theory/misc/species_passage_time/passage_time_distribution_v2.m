%% global settings
species_name = 'C3H8';
% end_t = '0.25718313951098054';
% end_t = '0.5';
% end_t = '0.9';
end_t = '1.0';

% text_str = 't=0.2 seconds';
% text_str = 't=0.5\tau';
% text_str = 't=0.9\tau';
text_str = 't=1.0\tau';

%% read file
file_dir_PT = fullfile(fileparts(mfilename('fullpath')));
% fn_PT = fullfile(file_dir_PT, 'data', 'species_pathway_AT_S62_HA3_0.25718313951098054.csv');
fn_PT = fullfile(file_dir_PT, 'data', ['species_pathway_AT_S62_HA3_', end_t ,'.csv']);

delimiter = '';
formatSpec = '%f%[^\n\r]';
%% Open the text file.
fileID = fopen(fn_PT,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN, 'ReturnOnError', false);
%% Close the text file.
fclose(fileID);
PT_vec = dataArray{:, 1};
%% Clear temporary variables
clearvars fn_PT delimiter formatSpec fileID dataArray ans;

%% passage time calculated from drc
fn_time_pt = fullfile(file_dir_PT, 'data', ['time_pt_dlsode_M_', end_t ,'.csv']);
delimiter = {''};
formatSpec = '%f%[^\n\r]';
fileID = fopen(fn_time_pt,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN,  'ReturnOnError', false);
fclose(fileID);
time_pt_vec = dataArray{:, 1};
clearvars fn_time_pt delimiter formatSpec fileID dataArray ans;

fn_drc_pt = fullfile(file_dir_PT, 'data', ['pt_dlsode_M_', end_t ,'.csv']);
delimiter = {''};
formatSpec = '%f%[^\n\r]';
fileID = fopen(fn_drc_pt,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN,  'ReturnOnError', false);
fclose(fileID);
pt_drc_vec = dataArray{:, 1};
clearvars fn_drc_pt delimiter formatSpec fileID dataArray ans;

%% plot
fig = figure();
% https://www.mathworks.com/help/matlab/graphics_transition/why-are-plot-lines-different-colors.html
% https://www.mathworks.com/help/matlab/creating_plots/customize-graph-with-two-y-axes.html
co = [    0    0.4470    0.7410 % 1th plot
%     0.8500    0.3250    0.0980 % 2nd plot
    0.9290    0.6940    0.1250 % 3rd plot
%     0.4940    0.1840    0.5560 % 4th plot
%     0.4660    0.6740    0.1880 % 5th plot
%     0.3010    0.7450    0.9330 % 6th plot
%     0.6350    0.0780    0.1840 % 7th plot
%     0         0    1.0000 % 8th plot, blue
    1   0   0 % for placeholder
    1   0   0]; % 9th plot, red
set(fig,'defaultAxesColorOrder',co)

% histogram
n_bins = 100;
p1 = histogram(PT_vec, n_bins, 'Normalization','pdf'); hold on;
% [N,edges] = histcounts(PT_vec, n_bins);
p2= plot(time_pt_vec, pt_drc_vec, 'color',co(2, :), 'LineWidth', 2); hold on;

%% settings
set(gca,'GridLineStyle','--');
xlabel('Time (seconds)', 'FontSize', 20);
ylabel('PDF ({second}^{-1})', 'FontSize', 20);
% ylim([0.0, 10^(0.9995*(log10(const_c1)))]);

leg_h = legend([p1; p2],... 
    'SOHR', 'EXACT');
set(leg_h, 'FontSize', 13, 'Box', 'off');
set(leg_h, 'Location', 'West')

%% global figure settings
grid on;
% xlim([1, end_idx1]);
% leg_h = legend([p1; p2; p3],'EXACT','SOHR', 'Percentage Error');
% set(leg_h, 'FontSize', 14, 'Box', 'off');
% % [left, bottom, weight, height]
% set(leg_h, 'Position', [0.275, 0.19, 0.5, 0.1])
% 
a_x = gca;
t_x = a_x.XLim(1) + 0.325*(a_x.XLim(2) - a_x.XLim(1));
t_y = a_x.YLim(1) + 0.818*(a_x.YLim(2) - a_x.YLim(1));
text(t_x, t_y, text_str, 'FontSize', 20);

%% Zoom in figure1
% % create a new pair of axes inside current figure
% z1_position = [.23 .45 .25 .25];
% axes('position',z1_position);
% box on; % put box around new pair of axes
% plot([path_idx(1), path_idx(end_idx2)], [const_c1, const_c1], ...
%     'color',co(1, :), 'LineWidth', 2); hold on;
% % cumulative
% semilogy(path_idx(1:delta2:end_idx2), cumu_path_p_vec1(1:delta2:end_idx2) , ...
%     'color',co(2, :), 'linestyle', '-', 'LineWidth', 2); hold on;
% axis tight;
% grid on;
% ylim([10^(1.0*log10(cumu_path_p_vec1(1))), 10^(0.999*(log10(const_c1)))]);


%% save to file
figname = strcat(species_name, '_species_passage_time_', end_t, '_v2.png');
print(fig, fullfile(file_dir_PT, figname), '-r200', '-dpng');
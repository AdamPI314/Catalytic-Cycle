%% global settings
spe_idx=17;
spe_name='CH_2O';
folder1 = '0.5tau';
tau = '0.777660157519';
end_t = '0.5';

%% Current file directory
file_dir = fullfile(fileparts(mfilename('fullpath')));

%% const concentration at a time
fn_conc_c1 = fullfile(file_dir, folder1, strcat('const_conc_500_',num2str(spe_idx),'_',tau,'_',end_t,'.csv'));
delimiter = '';
formatSpec = '%f%[^\n\r]';
%% Open the text file.
fileID = fopen(fn_conc_c1,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);
%% Close the text file.
fclose(fileID);
%% Allocate imported array to column variable names
const_c1 = dataArray{:, 1};
%% Clear temporary variables
clearvars fn_conc_c1 delimiter formatSpec fileID dataArray ans;

%% sorted pathway probabilities
file_dir = fullfile(fileparts(mfilename('fullpath')));
fn_path_p1 = fullfile(file_dir, folder1, strcat('path_prob_top_n_sorted_500_',num2str(spe_idx),'_',tau,'_',end_t,'.csv'));

delimiter = '';
formatSpec = '%f%[^\n\r]';
%% Open the text file.
fileID = fopen(fn_path_p1,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);
%% Close the text file.
fclose(fileID);
path_p_vec1 = dataArray{:, 1};
%% Clear temporary variables
clearvars fn_path_p1 delimiter formatSpec fileID dataArray ans;

%% cumulative pathway probabilities
cumu_path_p_vec1 = path_p_vec1;
for i = 2:length(cumu_path_p_vec1)
    cumu_path_p_vec1(i) = cumu_path_p_vec1(i) + cumu_path_p_vec1(i-1);
end

%% linear array
path_idx = linspace(1, length(path_p_vec1), length(path_p_vec1));
path_idx = path_idx';

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


% plot pathway probabilities
% individual
end_idx = 500;
delta_1 = 1;
p1 = semilogy(path_idx(1:delta_1:end_idx), path_p_vec1(1:delta_1:end_idx), ...
    'color',co(1, :), 'LineWidth', 2); hold on;
% cumulative
p2 = semilogy(path_idx(1:delta_1:end_idx), cumu_path_p_vec1(1:delta_1:end_idx) , ...
    'color',co(2, :), 'linestyle', '-', 'LineWidth', 2); hold on;
% placeholder for legend
p3 = plot(NaN, NaN, 'LineStyle', ':', 'LineWidth', 2, 'color', co(3, :));
%% conc
set(gca,'GridLineStyle','--');
xlabel('#Path', 'FontSize', 20);
ylabel('Probability', 'FontSize', 20);
% ylim([10^-8.3, 10^-2.7]);
ylim([10^(1.0*log10(path_p_vec1(end))), 10^(0.98*(log10(const_c1)))]);

yyaxis right
error_data1 = (const_c1 - cumu_path_p_vec1)./const_c1.*100;
semilogy(path_idx(1:delta_1:end_idx), error_data1(1:delta_1:end_idx), ... 
    'LineStyle', ':', 'LineWidth', 2, 'color', co(3, :)); hold on;
ylabel('%Relative Error', 'FontSize', 20);
% set(gca, 'ytick', []);

%% global settings
grid on;
xlim([1, end_idx]);
leg_h = legend([p1; p2; p3],'Pathway Probability','Cumulative Pathway Probability', 'Percentage Error');
set(leg_h, 'FontSize', 12, 'Box', 'off');
% [left, bottom, weight, height]
set(leg_h, 'Position', [0.4, 0.32, 0.5, 0.1])

%% Zoom in figure1
% create a new pair of axes inside current figure
axes('position',[.2 .5 .26 .35])
box on % put box around new pair of axes
delta_2 = 1;
end_idx2 = 10;
semilogy(path_idx(1:delta_2:end_idx2), path_p_vec1(1:delta_2:end_idx2), ...
    'marker', 'o', 'color',co(1, :), 'LineWidth', 2); hold on;
% cumulative
semilogy(path_idx(1:delta_2:end_idx2), cumu_path_p_vec1(1:delta_2:end_idx2) , ...
    'marker', 'o', 'color',co(2, :), 'linestyle', '-', 'LineWidth', 2); hold on;
axis tight;
grid on;
%% Error
yyaxis right
semilogy(path_idx(1:delta_2:end_idx2), error_data1(1:delta_2:end_idx2), ... 
    'marker', 'o', 'LineStyle', ':', 'LineWidth', 2, 'color', 'r'); hold on;
a_x = gca;
t_x = a_x.XLim(1) + 0.465*(a_x.XLim(2) - a_x.XLim(1));
t_y = a_x.YLim(1) + 0.618*(a_x.YLim(2) - a_x.YLim(1));
text(t_x, t_y, spe_name);

%% Zoom in figure2
co = [    0    0.4470    0.7410 % 1th plot
%     0.8500    0.3250    0.0980 % 2nd plot
%     0.9290    0.6940    0.1250 % 3rd plot
%     0.4940    0.1840    0.5560 % 4th plot
%     0.4660    0.6740    0.1880 % 5th plot
%     0.3010    0.7450    0.9330 % 6th plot
%     0.6350    0.0780    0.1840 % 7th plot
%     0         0    1.0000 % 8th plot, blue
    1   0   0 % for placeholder
    1   0   0]; % 9th plot, red
set(fig,'defaultAxesColorOrder',co);
z2_position = [.55 .5 .26 .35];
% create a new pair of axes inside current figure
axes('position',z2_position);

% plot a placeholder
plot(NaN, NaN);

set(gca, 'ytick', []);
yyaxis right
box on; % put box around new pair of axes

delta_3 = 1;
end_idx3 = 100;
zf2_1 = plot([path_idx(1), path_idx(end_idx3)], [const_c1, const_c1], ...
    'marker', 'o', 'color',co(1, :), 'LineWidth', 2); hold on;
% cumulative
zf2_2 = plot(path_idx(1:delta_3:end_idx3), cumu_path_p_vec1(1:delta_3:end_idx3) , ...
    'marker', 'o', 'color',co(2, :), 'linestyle', '-', 'LineWidth', 2); hold on;
axis tight;
grid on;
% legend
leg_z2 = legend([zf2_1; zf2_2],'EXACT','PATHWAY');
set(leg_z2, 'FontSize', 10, 'Box', 'off');
% [left, bottom, weight, height]
set(leg_z2, 'Position', [z2_position(1) + z2_position(3)*0.4, ...
    z2_position(2) + z2_position(4)*0.1, 0.1, 0.1])

%% save to file
figname = strcat('pathway_prob_concentration_',num2str(spe_idx), '.png');
print(fig, fullfile(file_dir, folder1, figname), '-r200', '-dpng');
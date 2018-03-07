%% global settings
spe_idx=14;
tau = 0.777660157519;
end_t = 0.9;

%% Current file directory
file_dir = fullfile(fileparts(mfilename('fullpath')));

%% const concentration at a time
fn_conc_c1 = fullfile(file_dir, '0.5tau', strcat('const_conc_500_',num2str(14),'_0.777660157519.csv'));
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
fn_path_p1 = fullfile(file_dir, '0.5tau', strcat('path_prob_top_n_sorted_500_',num2str(14),'_0.777660157519.csv'));

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

%% plot
%% plot
fig = figure();
% https://www.mathworks.com/help/matlab/graphics_transition/why-are-plot-lines-different-colors.html
% https://www.mathworks.com/help/matlab/creating_plots/customize-graph-with-two-y-axes.html
co = [    0    0.4470    0.7410 % 1th plot
    0.8500    0.3250    0.0980 % 2nd plot
    0.9290    0.6940    0.1250 % 3rd plot
%     0.4940    0.1840    0.5560 % 4th plot
%     0.4660    0.6740    0.1880 % 5th plot
%     0.3010    0.7450    0.9330 % 6th plot
%     0.6350    0.0780    0.1840 % 7th plot
%     0         0    1.0000 % 8th plot, blue
    1   0   0]; % 9th plot, red
set(fig,'defaultAxesColorOrder',co)


% plot conc
% npropyloo/npropyl
p1 = semilogy(time_vec, conc_mat(:, 79)./conc_mat(:, 61), 'color',co(1, :), 'LineWidth', 2); hold on;
% R_npropyloo * [npropyl] / (R_npropyl* [npropyloo])
p2 = semilogy(time_vec, (reaction_R_mat(:, 1069+1).*conc_mat(:, 79)) ./ (reaction_R_mat(:, 1068+1).*conc_mat(:, 61)) , ...
    'color',co(1, :), 'linestyle', '-.', 'LineWidth', 2); hold on;
% npropyloo/QOOH_1
p3 = semilogy(time_vec, conc_mat(:, 79)./conc_mat(:, 88), 'color',co(2, :), 'LineWidth', 2); hold on;

%% conc
set(gca,'GridLineStyle','--');
xlabel('Time (seconds)', 'FontSize', 20);
ylabel('Ratios', 'FontSize', 20);
ylim([10^1.7, 10^4.95]);

%% global settings
grid on;
xlim([0, tau*end_t]);
leg_h = legend([p1; p2; p3],'[nROO]/[nR]','K_{eq}','[nROO]/[QOOH_1]');
set(leg_h, 'FontSize', 14, 'Box', 'off');
set(leg_h, 'Location', 'South')


%% save to file
figname = strcat('pathway_prob_concentration_',num2str(spe_idx), '.png');
print(fig, fullfile(file_dir, 'output', figname), '-r200', '-dpng');
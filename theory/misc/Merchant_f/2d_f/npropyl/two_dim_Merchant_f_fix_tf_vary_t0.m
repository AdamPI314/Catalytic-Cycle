%% global settings
file_dir = fullfile(fileparts(mfilename('fullpath')));

spe_idx = '60';
atom_f = 'HA6';
spe_name = 'npropyl';
tau = 0.777660157519;
end_t = '0.9';
% end_t = '0.12859156975';
cycle = 'primary_cycle';
% cycle = 'others';
n_path = 100;

fn_2d_f = fullfile(file_dir, ['Merchant_f_2d_S', spe_idx, '_', atom_f, '_', end_t ,'.csv']);

delimiter = ',';
formatStr = "%f%f%f";
for i=1:n_path
    formatStr = formatStr + "%f";
end
formatStr = formatStr + "%[^\n\r]";
formatSpec = char(formatStr);

%% Open the text file.
fileID = fopen(fn_2d_f,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN,  'ReturnOnError', false);
%% Close the text file.
fclose(fileID);
f_mat = [dataArray{1:end-1}];

t0 = f_mat(:, 1);
tf = f_mat(:, 2);
% f_value = f_mat(:, end);

for i = 1:length(t0)
    t0(i) = t0(i) * tau;
    tf(i) = tf(i) * tau;
end

% path index
offset = 2;
path_idx1 = [1, 2, 3];
% path_idx1 = [4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100];
for i = 1:length(path_idx1)
    if i==1
        f_value1 = f_mat(:, offset + path_idx1(i));
    else
        f_value1 = f_value1 + f_mat(:, offset + path_idx1(i));
    end    
end
% construct 3d surface
xlin1 = linspace(min(t0), max(t0), 25);
ylin1 = linspace(min(tf), max(tf), 25);
[X1,Y1] = meshgrid(xlin1, ylin1);
f1 = scatteredInterpolant(t0, tf, f_value1);
% fix t_0, change t_f
X_tmp1 = X1(1, :);
Y_tmp1 = Y1(end, :);
Z_tmp1 = f1(X_tmp1, Y_tmp1);

% path_idx2 = [1, 2, 3];
path_idx2 = [4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100];
for i = 1:length(path_idx2)
    if i==1
        f_value2 = f_mat(:, offset + path_idx2(i));
    else
        f_value2 = f_value2 + f_mat(:, offset + path_idx2(i));
    end    
end
% construct 3d surface
xlin2 = linspace(min(t0), max(t0), 25);
ylin2 = linspace(min(tf), max(tf), 25);
[X2,Y2] = meshgrid(xlin2, ylin2);
f2 = scatteredInterpolant(t0, tf, f_value2);

% fix t_0, change t_f
X_tmp2 = X2(1, :);
Y_tmp2 = Y2(end, :);
Z_tmp2 = f2(X_tmp2, Y_tmp2);

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
    1   0   0 % placeholder
    1   0   0 % bl
    ]; 
set(fig,'defaultAxesColorOrder',co)

%% plot
p1 = plot(X_tmp1, Z_tmp1, 'LineWidth', 2, 'marker', '>'); hold on;
p2 = plot(X_tmp2, Z_tmp2, 'LineWidth', 2, 'marker', '<'); hold on;

%% settings
set(gca,'GridLineStyle','--');
xlabel('$t$ (seconds)', 'Interpreter','latex', 'FontSize', 20);
ylabel('$f$', 'Interpreter','latex', 'FontSize', 20);
% ylim([0.15, 1.8]);

%% global settings
grid on;
leg_h = legend([p1, p2], 'primary cycle', 'others');
set(leg_h, 'FontSize', 14, 'Box', 'off');
leg_pos = [0.5 0.325 0.1 0.1];
set(leg_h, 'Position', [leg_pos(1),leg_pos(2), leg_pos(3), leg_pos(4)]);

%% text
a_x = gca;
t_x = a_x.XLim(1) + 0.325*(a_x.XLim(2) - a_x.XLim(1));
t_y = a_x.YLim(1) + 0.808*(a_x.YLim(2) - a_x.YLim(1));
text(t_x, t_y, 'npropyl, $t_f$=0.9$\tau$', 'Interpreter','latex', 'FontSize', 20);

%% save to file
figname = strcat('1d_Merchant_f_fix_tf_vary_t0_', end_t, '_S', spe_idx, '_', cycle, '_v5.png');
print(fig, fullfile(file_dir, figname), '-r200', '-dpng');

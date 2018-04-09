%% global settings
file_dir = fullfile(fileparts(mfilename('fullpath')));

spe_idx = '60';
atom_f = 'HA6';
spe_name = 'npropyl';
tau = 0.777660157519;
end_t = '0.9';
% end_t = '0.12859156975';
cycle = 'all';
% cycle = 'primary_cycle';
% cycle = 'others';
n_path = 1000;

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
path_idx = linspace(1, n_path, n_path-1+1);
% path_idx = [1, 2, 3];
% path_idx = linspace(3, n_path, n_path-3+1);
for i = 1:length(path_idx)
    if i==1
        f_value = f_mat(:, offset + path_idx(i));
    else
        f_value = f_value + f_mat(:, offset + path_idx(i));
    end    
end

% construct 3d surface
xlin = linspace(min(t0), max(t0), 25);
ylin = linspace(min(tf), max(tf), 25);
[X,Y] = meshgrid(xlin, ylin);
f = scatteredInterpolant(t0, tf, f_value);
Z = f(X,Y);

%% update Z
for i = 1:length(X)    
    for j = i+1:length(X)
        Z(i,j) = nan;
    end
end

%% plot
fig = figure();

%% plot
% mesh(X,Y,Z); %interpolated

% reduce number of s.f. in coutour plot
contour(X,Y,Z, 20, 'ShowText', 'on');
hold on;

% % draw circle to emphase
% he1 = ellipse(0.045,0.075,pi/2,0.085,0.65);
% set(he1, 'LineStyle', '--', 'color', 'b', 'LineWidth', 2.0);
% hold on;

% draw circle to emphase
% he2 = ellipse(0.065,0.075,pi/2,0.085,0.525);
% set(he2, 'LineStyle', '--', 'color', 'r', 'LineWidth', 2.0);
% hold on;

% center_x2 = 0.085; center_y2 = 0.525; delta_x2 = 0.065*2; delta_y2 = 0.075*2;
% r2=rectangle('Position',[center_x2 - delta_x2/2, center_y2 - delta_y2/2, delta_x2, delta_y2]);
% r2.EdgeColor = 'r';
% r2.LineWidth = 3;
% r2.LineStyle = '--';
% r2.Curvature = [0.5,1];
% hold on;

% % draw circle to emphase
% he3 = ellipse(0.11,0.065,pi/2,0.075,0.325);
% set(he3, 'LineStyle', '--', 'color', 'k', 'LineWidth', 2.0);
% hold on;

axis tight;
hold on;


%% settings
set(gca,'GridLineStyle','--');
xlabel('$t$ (seconds)', 'Interpreter','latex', 'FontSize', 20);
ylabel('$t_f$ (seconds)', 'Interpreter','latex', 'FontSize', 20);
% xlim([0, tau*str2double(end_t)]);
% ylim([0, tau*str2double(end_t)]);
grid on;

%% text
a_x = gca;
t_x = a_x.XLim(1) + 0.525*(a_x.XLim(2) - a_x.XLim(1));
t_y = a_x.YLim(1) + 0.278*(a_x.YLim(2) - a_x.YLim(1));
text(t_x, t_y, [spe_name, '@ $t$' char(10) 'stop path@ $t_f$' char(10) '$t_f > t$'], 'Interpreter','latex', 'FontSize', 20);

%% save to file
figname = strcat('2d_Merchant_f_', end_t, '_S', spe_idx, '_', cycle, '_v4.png');
print(fig, fullfile(file_dir, figname), '-r200', '-dpng');



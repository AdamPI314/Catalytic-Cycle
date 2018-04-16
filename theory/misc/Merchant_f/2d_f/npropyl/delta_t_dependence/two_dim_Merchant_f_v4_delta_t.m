%% global settings
file_dir = fullfile(fileparts(mfilename('fullpath')));

% marker
markers = {'+' , 'o' , '*' , 'x' , 'square' , 'diamond' , 'v' , '^' , '>' , '<' , 'pentagram' , 'hexagram' , '.', 'none'};
spe_idx = '60';
atom_f = 'HA6';
spe_name = 'npropyl';
tau = 0.777660157519;
end_t = '0.9';
% end_t = '0.12859156975';
% cycle = 'all';
% cycle = 'p1';
% cycle = 'p2';
% cycle = 'p3';
cycle = 'primary_cycle';
% cycle = 'others';
n_path = 3;

% fn_2d_f = fullfile(file_dir, ['Merchant_f_2d_S', spe_idx, '_', atom_f, '_', end_t ,'.csv']);
fn_2d_f = fullfile(file_dir, ['Merchant_f_2d', '.csv']);

delimiter = ',';
formatStr = '%f%f%f';
for i=1:n_path
    formatStr = strcat(formatStr, '%f');
end
formatStr = strcat(formatStr, '%[^\n\r]');
formatSpec = char(formatStr);

%% Open the text file.
fileID = fopen(fn_2d_f,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue', NaN,  'ReturnOnError', false);
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

tf = tf - t0;

% path index
offset = 2;
% path_idx = linspace(1, n_path, n_path-1+1);
% path_idx = [1];
% path_idx = [2];
% path_idx = [3];
path_idx = [1, 2, 3];
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
    for j = length(X) - i + 1 : length(X)
        Z(i,j) = nan;
    end
end

%% plot
fig = figure();

%% plot
X_tmp1 = X(1, :);

delta_t_vec = [1e-5, 5e-5, 1e-4, 5e-4, 1e-3, 2.5e-3, 5e-3, 1e-2, 1.5e-2, 2.5e-2, 5e-2, 1e-1];
N = length(delta_t_vec);
colors = lines(N);
% colors = colorcube(N);
str_name = cell(N,1);
for i=1:N
    str_name{i, 1} = strcat('$\delta$t=', num2str(delta_t_vec(i),'%1.1e\n'));
end

H = gobjects(N);

for idx=1:N
    % delta
    delta_t = ones(1, length(X_tmp1));
    delta_t = delta_t.* delta_t_vec(idx);
    Z_tmp1 = f(X_tmp1, delta_t);
    % check data
    for i=1:length(X_tmp1)
%         if X_tmp1(i) + delta_t(i) > ylin(end)
        if X_tmp1(i) + delta_t(i) > str2double(end_t)
            Z_tmp1(i) = nan;
        end
    end
    H(idx) = plot(X_tmp1, Z_tmp1, 'LineWidth', 2, 'color', colors(idx, :), ...
        'marker', markers{1, mod(idx-1, length(markers))+ 1}); hold on;
    hold on;
end




%% settings
set(gca,'GridLineStyle','--');
xlabel('$t$ (seconds)', 'Interpreter','latex', 'FontSize', 20);
ylabel('$\#$(OH)', 'Interpreter','latex', 'FontSize', 20);
% xlim([0, tau*str2double(end_t)]);
% ylim([0, tau*str2double(end_t)]);
grid on;

%%  legend
leg_h = legend(str_name, 'Interpreter','latex');
set(leg_h, 'FontSize', 12, 'Box', 'off');
% set(leg_h, 'Location', 'South');

% %% text
% a_x = gca;
% t_x = a_x.XLim(1) + 0.525*(a_x.XLim(2) - a_x.XLim(1));
% t_y = a_x.YLim(1) + 0.278*(a_x.YLim(2) - a_x.YLim(1));
% text(t_x, t_y, [spe_name, '@ $t$' char(10) 'stop path@ $t_f$' char(10) '$t_f > t$'], 'Interpreter','latex', 'FontSize', 20);

%% save to file
figname = strcat('2d_Merchant_f_', end_t, '_S', spe_idx, '_', cycle, '_v4_delta_t.png');
print(fig, fullfile(file_dir, figname), '-r200', '-dpng');



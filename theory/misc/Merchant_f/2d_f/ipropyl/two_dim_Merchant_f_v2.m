%% global settings
file_dir = fullfile(fileparts(mfilename('fullpath')));

spe_idx = '60';
spe_name = 'ipropyl';
tau = 0.777660157519;
end_t = '0.9';

fn_2d_f = fullfile(file_dir, ['Merchant_f_2d' ,'.csv']);
delimiter = ',';
formatSpec = '%f%f%f%[^\n\r]';
fileID = fopen(fn_2d_f,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);
fclose(fileID);
t0 = dataArray{:, 1};
tf = dataArray{:, 2};
f_value = dataArray{:, 3};
clearvars fn_2d_f delimiter formatSpec fileID dataArray ans;

for i = 1:length(t0)
    t0(i) = t0(i) * tau;
    tf(i) = tf(i) * tau;
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
contour(X,Y,Z, 15, 'ShowText', 'on');

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
text(t_x, t_y, [spe_name, '@ $t$' char(10) 'stop path@ $t_f$'], 'Interpreter','latex', 'FontSize', 20);

%% save to file
figname = strcat('2d_Merchant_f_', end_t, '_S', spe_idx, '_v2.png');
print(fig, fullfile(file_dir, figname), '-r200', '-dpng');




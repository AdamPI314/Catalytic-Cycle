%% global settings
file_dir = fullfile(fileparts(mfilename('fullpath')));

spe_idx = '60';
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
mesh(X,Y,Z); %interpolated
axis tight;
hold on;
% plot3(t0,tf,f_value,'.','MarkerSize',15); %nonuniform

%% settings
set(gca,'GridLineStyle','--');
xlabel('$t_0$ (seconds)', 'Interpreter','latex', 'FontSize', 20);
ylabel('$t_f$ (seconds)', 'Interpreter','latex', 'FontSize', 20);
zlabel('$f$', 'Interpreter','latex', 'FontSize', 20);

%% viewpoint
az = -37.5;
el = 30+25;
view(az, el);

%% save to file
figname = strcat('2d_Merchant_f_', end_t, '_S', spe_idx, '_v1.png');
print(fig, fullfile(file_dir, figname), '-r200', '-dpng');




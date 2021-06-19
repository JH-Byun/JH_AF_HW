load 2_c_100.mat
e_x_100 = X - xhat;
ESS_100 = ESS;
time_spent_100 = time_spent;

load 2_c_500.mat
e_x_500 = X - xhat;
ESS_500 = ESS;
time_spent_500 = time_spent;

load 2_c_1000.mat
e_x_1000 = X - xhat;
ESS_1000 = ESS;
time_spent_1000 = time_spent;

MC = size(e_x_100,1);
simLen = size(e_x_100,2);

for i = 1:simLen
    RMS_100(:,i) = sqrt(sum(e_x_100(:,i).^2)/MC);
    RMS_500(:,i) = sqrt(sum(e_x_500(:,i).^2)/MC);
    RMS_1000(:,i) = sqrt(sum(e_x_1000(:,i).^2)/MC);
end

% settings
default_font_size = 12.5; % default font size
default_line_width = 2; % default line width

% N = 100
figure('DefaultAxesFontSize',default_font_size);
subplot(2,1,1)
plot(1:simLen,RMS_100,'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
title('SIS: RMSE (N = 100)')
subplot(2,1,2)
plot(1:simLen,ESS_100,'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
title('SIS: ESS (N = 100)')

% UKF
figure('DefaultAxesFontSize',default_font_size);
subplot(2,1,1)
plot(1:simLen,RMS_500,'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
title('SIS: RMSE (N = 500)')
subplot(2,1,2)
plot(1:simLen,ESS_500,'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
title('SIS: ESS (N = 500)')

% SIS
figure('DefaultAxesFontSize',default_font_size);
subplot(2,1,1)
plot(1:simLen,RMS_1000,'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
title('SIS: RMSE (N = 1000)')
subplot(2,1,2)
plot(1:simLen,ESS_1000,'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
title('SIS: ESS (N = 1000)')

% computation time comparison
figure('DefaultAxesFontSize',default_font_size);
plot(1:MC, time_spent_100,'LineWidth',default_line_width)
hold on
plot(1:MC, time_spent_500,'LineWidth',default_line_width)
hold on
plot(1:MC, time_spent_1000,'LineWidth',default_line_width)
grid on
xlabel('MC iteration step [-]')
title('Computation time comparison: particle numbers N')
legend({'N = 100','N = 500','N = 1000'},'Location', 'best')
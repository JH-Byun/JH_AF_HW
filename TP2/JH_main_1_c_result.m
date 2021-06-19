load 1_c_100.mat
e_x_100 = x_true - xhat_SIS;
sgmhat_100 = sgmhat_SIS;
time_spent_100 = time_spent_SIS;

load 1_c_1000.mat
e_x_1000 = x_true - xhat_SIS;
sgmhat_1000 = sgmhat_SIS;
time_spent_1000 = time_spent_SIS;

load 1_c_10000.mat
e_x_10000 = x_true - xhat_SIS;
sgmhat_10000 = sgmhat_SIS;
time_spent_10000 = time_spent_SIS;

MC = size(e_x_100,1)/2;
simLen = size(e_x_100,2);

for i = 1:simLen
    RMS_100(:,i) = [sqrt(sum(e_x_100(1:2:end,i).^2)/MC);sqrt(sum(e_x_100(2:2:end,i).^2)/MC)];
    RMS_1000(:,i) = [sqrt(sum(e_x_1000(1:2:end,i).^2)/MC);sqrt(sum(e_x_1000(2:2:end,i).^2)/MC)];
    RMS_10000(:,i) = [sqrt(sum(e_x_10000(1:2:end,i).^2)/MC);sqrt(sum(e_x_10000(2:2:end,i).^2)/MC)]; 
end

% settings
default_font_size = 12.5; % default font size
default_line_width = 2; % default line width

% case 1
figure('DefaultAxesFontSize',default_font_size);
subplot(2,1,1)
plot(1:simLen,RMS_100(1,:),'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
title('Position [m]')
subplot(2,1,2)
plot(1:simLen,RMS_100(2,:),'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
title('Velocity [m/s]')
sgtitle('SIS RMSE result (N = 100)')

figure;
subplot(2,1,1)
plot(1:simLen,sgmhat_100(1:2:end,:),'LineWidth',default_line_width)
grid on
subplot(2,1,2)
plot(1:simLen,sgmhat_100(2:2:end,:),'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
sgtitle('SIS root(P) result (N = 100)')

% case 2
figure('DefaultAxesFontSize',default_font_size);
subplot(2,1,1)
plot(1:simLen,RMS_1000(1,:),'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
title('Position [m]')
subplot(2,1,2)
plot(1:simLen,RMS_1000(2,:),'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
title('Velocity [m/s]')
sgtitle('SIS RMSE result (N = 1000)')

figure;
subplot(2,1,1)
plot(1:simLen,sgmhat_1000(1:2:end,:),'LineWidth',default_line_width)
grid on
subplot(2,1,2)
plot(1:simLen,sgmhat_1000(2:2:end,:),'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
sgtitle('SIS root(P) result (N = 1000)')

% case 3
figure('DefaultAxesFontSize',default_font_size);
subplot(2,1,1)
plot(1:simLen,RMS_10000(1,:),'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
title('Position [m]')
subplot(2,1,2)
plot(1:simLen,RMS_10000(2,:),'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
title('Velocity [m/s]')
sgtitle('SIS RMSE result (N = 10000)')

figure;
subplot(2,1,1)
plot(1:simLen,sgmhat_10000(1:2:end,:),'LineWidth',default_line_width)
grid on
subplot(2,1,2)
plot(1:simLen,sgmhat_10000(2:2:end,:),'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
sgtitle('SIS root(P) result (N = 10000)')

%% computation time
figure('DefaultAxesFontSize',default_font_size);
plot(1:MC,time_spent_100,'LineWidth',default_line_width)
hold on
plot(1:MC,time_spent_1000,'LineWidth',default_line_width)
hold on
plot(1:MC,time_spent_10000,'LineWidth',default_line_width)
grid on
xlabel('Monte-Carlo iteration step [-]')
legend({'N = 100','N = 1000','N = 10000'},'Location','best')
title('SIS computation result')

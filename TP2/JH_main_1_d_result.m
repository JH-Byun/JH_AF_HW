load 1_d_1.mat
e_x_SIS_1 = x_true - xhat_SIS;
sgmhat_SIS_1 = sgmhat_SIS;
time_spent_SIS_1 = time_spent_SIS;
e_x_SIR_1 = x_true - xhat_SIR;
sgmhat_SIR_1 = sgmhat_SIR;
time_spent_SIR_1 = time_spent_SIR;

load 1_d_2.mat
e_x_SIS_2 = x_true - xhat_SIS;
sgmhat_SIS_2 = sgmhat_SIS;
time_spent_SIS_2 = time_spent_SIS;
e_x_SIR_2 = x_true - xhat_SIR;
sgmhat_SIR_2 = sgmhat_SIR;
time_spent_SIR_2 = time_spent_SIR;

load 1_d_3.mat
e_x_SIS_3 = x_true - xhat_SIS;
sgmhat_SIS_3 = sgmhat_SIS;
time_spent_SIS_3 = time_spent_SIS;
e_x_SIR_3 = x_true - xhat_SIR;
sgmhat_SIR_3 = sgmhat_SIR;
time_spent_SIR_3 = time_spent_SIR;

load 1_d_4.mat
e_x_SIS_4 = x_true - xhat_SIS;
sgmhat_SIS_4 = sgmhat_SIS;
time_spent_SIS_4 = time_spent_SIS;
e_x_SIR_4 = x_true - xhat_SIR;
sgmhat_SIR_4 = sgmhat_SIR;
time_spent_SIR_4 = time_spent_SIR;

MC = size(e_x_SIS_1,1)/2;
simLen = size(e_x_SIS_1,2);

for i = 1:simLen
    RMS_SIS_1(:,i) = [sqrt(sum(e_x_SIS_1(1:2:end,i).^2)/MC);sqrt(sum(e_x_SIS_1(2:2:end,i).^2)/MC)];
    RMS_SIS_2(:,i) = [sqrt(sum(e_x_SIS_2(1:2:end,i).^2)/MC);sqrt(sum(e_x_SIS_2(2:2:end,i).^2)/MC)];
    RMS_SIS_3(:,i) = [sqrt(sum(e_x_SIS_3(1:2:end,i).^2)/MC);sqrt(sum(e_x_SIS_3(2:2:end,i).^2)/MC)];
    RMS_SIS_4(:,i) = [sqrt(sum(e_x_SIS_4(1:2:end,i).^2)/MC);sqrt(sum(e_x_SIS_4(2:2:end,i).^2)/MC)];
    RMS_SIR_1(:,i) = [sqrt(sum(e_x_SIR_1(1:2:end,i).^2)/MC);sqrt(sum(e_x_SIR_1(2:2:end,i).^2)/MC)];
    RMS_SIR_2(:,i) = [sqrt(sum(e_x_SIR_2(1:2:end,i).^2)/MC);sqrt(sum(e_x_SIR_2(2:2:end,i).^2)/MC)];
    RMS_SIR_3(:,i) = [sqrt(sum(e_x_SIR_3(1:2:end,i).^2)/MC);sqrt(sum(e_x_SIR_3(2:2:end,i).^2)/MC)];
    RMS_SIR_4(:,i) = [sqrt(sum(e_x_SIR_4(1:2:end,i).^2)/MC);sqrt(sum(e_x_SIR_4(2:2:end,i).^2)/MC)];
end

% settings
default_font_size = 12.5; % default font size
default_line_width = 2; % default line width

% case 1
figure('DefaultAxesFontSize',default_font_size);
subplot(2,1,1)
plot(1:simLen,RMS_SIS_1(1,:),'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
title('Position [m]')
subplot(2,1,2)
plot(1:simLen,RMS_SIS_1(2,:),'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
title('Velocity [m/s]')
sgtitle('SIS RMSE result (q = 10, r = 1)')

figure;
subplot(2,1,1)
plot(1:simLen,sgmhat_SIS_1(1:2:end,:),'LineWidth',default_line_width)
grid on
subplot(2,1,2)
plot(1:simLen,sgmhat_SIS_1(2:2:end,:),'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
sgtitle('SIS root(P) result (q = 10, r = 1)')

figure('DefaultAxesFontSize',default_font_size);
subplot(2,1,1)
plot(1:simLen,RMS_SIR_1(1,:),'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
title('Position [m]')
subplot(2,1,2)
plot(1:simLen,RMS_SIR_1(2,:),'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
title('Velocity [m/s]')
sgtitle('SIR RMSE result (q = 10, r = 1)')

figure;
subplot(2,1,1)
plot(1:simLen,sgmhat_SIR_1(1:2:end,:),'LineWidth',default_line_width)
grid on
subplot(2,1,2)
plot(1:simLen,sgmhat_SIR_1(2:2:end,:),'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
sgtitle('SIR root(P) result (q = 10, r = 1)')

% case 2
figure('DefaultAxesFontSize',default_font_size);
subplot(2,1,1)
plot(1:simLen,RMS_SIS_2(1,:),'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
title('Position [m]')
subplot(2,1,2)
plot(1:simLen,RMS_SIS_2(2,:),'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
title('Velocity [m/s]')
sgtitle('SIS RMSE result (q = 1, r = 1)')

figure;
subplot(2,1,1)
plot(1:simLen,sgmhat_SIS_2(1:2:end,:),'LineWidth',default_line_width)
grid on
subplot(2,1,2)
plot(1:simLen,sgmhat_SIS_2(2:2:end,:),'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
sgtitle('SIS root(P) result (q = 1, r = 1)')

figure('DefaultAxesFontSize',default_font_size);
subplot(2,1,1)
plot(1:simLen,RMS_SIR_2(1,:),'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
title('Position [m]')
subplot(2,1,2)
plot(1:simLen,RMS_SIR_2(2,:),'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
title('Velocity [m/s]')
sgtitle('SIR RMSE result (q = 1, r = 1)')

figure;
subplot(2,1,1)
plot(1:simLen,sgmhat_SIR_2(1:2:end,:),'LineWidth',default_line_width)
grid on
subplot(2,1,2)
plot(1:simLen,sgmhat_SIR_2(2:2:end,:),'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
sgtitle('SIR root(P) result (q = 1, r = 1)')

% case 3
figure('DefaultAxesFontSize',default_font_size);
subplot(2,1,1)
plot(1:simLen,RMS_SIS_3(1,:),'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
title('Position [m]')
subplot(2,1,2)
plot(1:simLen,RMS_SIS_3(2,:),'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
title('Velocity [m/s]')
sgtitle('SIS RMSE result (q = 0.1, r = 1)')

figure;
subplot(2,1,1)
plot(1:simLen,sgmhat_SIS_3(1:2:end,:),'LineWidth',default_line_width)
grid on
subplot(2,1,2)
plot(1:simLen,sgmhat_SIS_3(2:2:end,:),'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
sgtitle('SIS root(P) result (q = 0.1, r = 1)')

figure('DefaultAxesFontSize',default_font_size);
subplot(2,1,1)
plot(1:simLen,RMS_SIR_3(1,:),'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
title('Position [m]')
subplot(2,1,2)
plot(1:simLen,RMS_SIR_3(2,:),'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
title('Velocity [m/s]')
sgtitle('SIR RMSE result (q = 0.1, r = 1)')

figure;
subplot(2,1,1)
plot(1:simLen,sgmhat_SIR_3(1:2:end,:),'LineWidth',default_line_width)
grid on
subplot(2,1,2)
plot(1:simLen,sgmhat_SIR_3(2:2:end,:),'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
sgtitle('SIR root(P) result (q = 0.1, r = 1)')

% case 4
figure('DefaultAxesFontSize',default_font_size);
subplot(2,1,1)
plot(1:simLen,RMS_SIS_4(1,:),'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
title('Position [m]')
subplot(2,1,2)
plot(1:simLen,RMS_SIS_4(2,:),'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
title('Velocity [m/s]')
sgtitle('SIS RMSE result (q = 0.01, r = 1)')

figure;
subplot(2,1,1)
plot(1:simLen,sgmhat_SIS_4(1:2:end,:),'LineWidth',default_line_width)
grid on
subplot(2,1,2)
plot(1:simLen,sgmhat_SIS_4(2:2:end,:),'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
sgtitle('SIS root(P) result (q = 0.01, r = 1)')

figure('DefaultAxesFontSize',default_font_size);
subplot(2,1,1)
plot(1:simLen,RMS_SIR_4(1,:),'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
title('Position [m]')
subplot(2,1,2)
plot(1:simLen,RMS_SIR_4(2,:),'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
title('Velocity [m/s]')
sgtitle('SIR RMSE result (q = 0.01, r = 1)')

figure;
subplot(2,1,1)
plot(1:simLen,sgmhat_SIR_4(1:2:end,:),'LineWidth',default_line_width)
grid on
subplot(2,1,2)
plot(1:simLen,sgmhat_SIR_4(2:2:end,:),'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
sgtitle('SIR root(P) result (q = 0.01, r = 1)')

%% computation time
figure('DefaultAxesFontSize',default_font_size);
plot(1:MC,time_spent_SIS_1,'LineWidth',default_line_width)
hold on
plot(1:MC,time_spent_SIR_1,'LineWidth',default_line_width)
grid on
xlabel('Monte-Carlo iteration step [-]')
legend({'SIS','SIR'},'Location','best')
title('q = 10, r = 1: computation result')
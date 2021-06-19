load 2_a_EKF.mat
e_x_EKF = X - xhat;
sgmhat_EKF = sgmhat;
time_spent_EKF = time_spent;

load 2_a_UKF.mat
e_x_UKF = X - xhat;
sgmhat_UKF = sgmhat;
time_spent_UKF = time_spent;

load 2_a_SIS.mat
e_x_SIS = X - xhat;
sgmhat_SIS = sgmhat;
time_spent_SIS = time_spent;

load 2_a_SIR.mat
e_x_SIR = X - xhat;
sgmhat_SIR = sgmhat;
time_spent_SIR = time_spent;

MC = size(e_x_EKF,1);
simLen = size(e_x_EKF,2);

for i = 1:simLen
    RMS_EKF(:,i) = sqrt(sum(e_x_EKF(:,i).^2)/MC);
    RMS_UKF(:,i) = sqrt(sum(e_x_UKF(:,i).^2)/MC);
    RMS_SIS(:,i) = sqrt(sum(e_x_SIS(:,i).^2)/MC);
    RMS_SIR(:,i) = sqrt(sum(e_x_SIR(:,i).^2)/MC);
end

% settings
default_font_size = 12.5; % default font size
default_line_width = 2; % default line width

% EKF
figure('DefaultAxesFontSize',default_font_size);
subplot(2,1,1)
plot(1:simLen,RMS_EKF,'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
title('EKF: RMSE')
subplot(2,1,2)
plot(1:simLen,sgmhat_EKF,'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
title('EKF: root(P)')

% UKF
figure('DefaultAxesFontSize',default_font_size);
subplot(2,1,1)
plot(1:simLen,RMS_UKF,'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
title('UKF: RMSE')
subplot(2,1,2)
plot(1:simLen,sgmhat_UKF,'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
title('UKF: root(P)')

% SIS
figure('DefaultAxesFontSize',default_font_size);
subplot(2,1,1)
plot(1:simLen,RMS_SIS,'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
title('SIS: RMSE')
subplot(2,1,2)
plot(1:simLen,sgmhat_SIS,'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
title('SIS: root(P)')

% SIR
figure('DefaultAxesFontSize',default_font_size);
subplot(2,1,1)
plot(1:simLen,RMS_SIR,'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
title('SIR: RMSE')
subplot(2,1,2)
plot(1:simLen,sgmhat_SIR,'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
title('SIR: root(P)')

% comparison: EKF, UKF, SIS, SIR
figure('DefaultAxesFontSize',default_font_size);
plot(1:simLen, RMS_EKF,'LineWidth',default_line_width)
hold on
plot(1:simLen, RMS_UKF,'LineWidth',default_line_width)
hold on
plot(1:simLen, RMS_SIS,'LineWidth',default_line_width)
hold on
plot(1:simLen, RMS_SIR,'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
title('RMSE comparison: EKF, UKF, SIS, SIR')
legend({'EKF','UKF','SIS','SIR'},'Location', 'best')

% computation time comparison
figure('DefaultAxesFontSize',default_font_size);
plot(1:MC, time_spent_EKF,'LineWidth',default_line_width)
hold on
plot(1:MC, time_spent_UKF,'LineWidth',default_line_width)
hold on
plot(1:MC, time_spent_SIS,'LineWidth',default_line_width)
hold on
plot(1:MC, time_spent_SIR,'LineWidth',default_line_width)
grid on
xlabel('MC iteration step [-]')
title('Computation time comparison: EKF, UKF, SIS, SIR')
legend({'EKF','UKF','SIS','SIR'},'Location', 'best')
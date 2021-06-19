load 1_e_PF_SIS_1.mat
e_x_SIS_1 = x_true - xhat_SIS;
sgmhat_SIS_1 = sgmhat_SIS;
time_spent_SIS_1 = time_spent_SIS;
ESS_SIS_1 = ESS;

load 1_e_PF_SIS_2.mat
e_x_SIS_2 = x_true - xhat_SIS;
sgmhat_SIS_2 = sgmhat_SIS;
time_spent_SIS_2 = time_spent_SIS;
ESS_SIS_2 = ESS;

load 1_e_PF_SIS_3.mat
e_x_SIS_3 = x_true - xhat_SIS;
sgmhat_SIS_3 = sgmhat_SIS;
time_spent_SIS_3 = time_spent_SIS;
ESS_SIS_3 = ESS;

load 1_e_PF_SIS_4.mat
e_x_SIS_4 = x_true - xhat_SIS;
sgmhat_SIS_4 = sgmhat_SIS;
time_spent_SIS_4 = time_spent_SIS;
ESS_SIS_4 = ESS;

MC = size(e_x_SIS_1,1)/2;
simLen = size(e_x_SIS_1,2);

for i = 1:simLen
    RMS_SIS_1(:,i) = [sqrt(sum(e_x_SIS_1(1:2:end,i).^2)/MC);sqrt(sum(e_x_SIS_1(2:2:end,i).^2)/MC)];
    RMS_SIS_2(:,i) = [sqrt(sum(e_x_SIS_2(1:2:end,i).^2)/MC);sqrt(sum(e_x_SIS_2(2:2:end,i).^2)/MC)];
    RMS_SIS_3(:,i) = [sqrt(sum(e_x_SIS_3(1:2:end,i).^2)/MC);sqrt(sum(e_x_SIS_3(2:2:end,i).^2)/MC)];
    RMS_SIS_4(:,i) = [sqrt(sum(e_x_SIS_4(1:2:end,i).^2)/MC);sqrt(sum(e_x_SIS_4(2:2:end,i).^2)/MC)];
end

% settings
default_font_size = 12.5; % default font size
default_line_width = 2; % default line width

% case 1
figure('DefaultAxesFontSize',default_font_size);
subplot(2,1,1)
plot(1:simLen,RMS_SIS_1(1,:),'LineWidth',default_line_width)
hold on
plot(1:simLen,RMS_SIS_1(2,:),'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
legend({'RMSE of position','RMSE of velocity'},'Location', 'best')
subplot(2,1,2)
plot(1:simLen,ESS_SIS_1,'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
ylabel('ESS [-]')
sgtitle('SIS result (q = 10, r = 5)')

% case 2
figure('DefaultAxesFontSize',default_font_size);
subplot(2,1,1)
plot(1:simLen,RMS_SIS_2(1,:),'LineWidth',default_line_width)
hold on
plot(1:simLen,RMS_SIS_2(2,:),'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
legend({'RMSE of position','RMSE of velocity'},'Location', 'best')
subplot(2,1,2)
plot(1:simLen,ESS_SIS_2,'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
ylabel('ESS [-]')
sgtitle('SIS result (q = 1, r = 5)')

% case 3
figure('DefaultAxesFontSize',default_font_size);
subplot(2,1,1)
plot(1:simLen,RMS_SIS_3(1,:),'LineWidth',default_line_width)
hold on
plot(1:simLen,RMS_SIS_3(2,:),'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
legend({'RMSE of position','RMSE of velocity'},'Location', 'best')
subplot(2,1,2)
plot(1:simLen,ESS_SIS_3,'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
ylabel('ESS [-]')
sgtitle('SIS result (q = 0.1, r = 5)')

% case 4
figure('DefaultAxesFontSize',default_font_size);
subplot(2,1,1)
plot(1:simLen,RMS_SIS_4(1,:),'LineWidth',default_line_width)
hold on
plot(1:simLen,RMS_SIS_4(2,:),'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
legend({'RMSE of position','RMSE of velocity'},'Location', 'best')
subplot(2,1,2)
plot(1:simLen,ESS_SIS_4,'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
ylabel('ESS [-]')
sgtitle('SIS result (q = 0.01, r = 5)')
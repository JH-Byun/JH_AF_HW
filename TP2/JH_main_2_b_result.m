load 2_b_1.mat
e_x_1 = X - xhat;
sgmhat_1 = sgmhat;
time_spent_1 = time_spent;

load 2_b_2.mat
e_x_2 = X - xhat;
sgmhat_2 = sgmhat;
time_spent_2 = time_spent;

load 2_b_3.mat
e_x_3 = X - xhat;
sgmhat_3 = sgmhat;
time_spent_3 = time_spent;

load 2_b_4.mat
e_x_4 = X - xhat;
sgmhat_4 = sgmhat;
time_spent_4 = time_spent;

load 2_b_5.mat
e_x_5 = X - xhat;
sgmhat_5 = sgmhat;
time_spent_5 = time_spent;

load 2_b_6.mat
e_x_6 = X - xhat;
sgmhat_6 = sgmhat;
time_spent_6 = time_spent;

load 2_b_7.mat
e_x_7 = X - xhat;
sgmhat_7 = sgmhat;
time_spent_7 = time_spent;

load 2_b_8.mat
e_x_8 = X - xhat;
sgmhat_8 = sgmhat;
time_spent_8 = time_spent;

load 2_b_9.mat
e_x_9 = X - xhat;
sgmhat_9 = sgmhat;
time_spent_9 = time_spent;

load 2_b_10.mat
e_x_10 = X - xhat;
sgmhat_10 = sgmhat;
time_spent_10 = time_spent;

load 2_b_11.mat
e_x_11 = X - xhat;
sgmhat_11 = sgmhat;
time_spent_11 = time_spent;

load 2_b_12.mat
e_x_12 = X - xhat;
sgmhat_12 = sgmhat;
time_spent_12 = time_spent;

MC = size(e_x_1,1);
simLen = size(e_x_1,2);

for i = 1:simLen
    RMS_1(:,i) = sqrt(sum(e_x_1(:,i).^2)/MC);
    RMS_2(:,i) = sqrt(sum(e_x_2(:,i).^2)/MC);
    RMS_3(:,i) = sqrt(sum(e_x_3(:,i).^2)/MC);
    RMS_4(:,i) = sqrt(sum(e_x_4(:,i).^2)/MC);
    RMS_5(:,i) = sqrt(sum(e_x_5(:,i).^2)/MC);
    RMS_6(:,i) = sqrt(sum(e_x_6(:,i).^2)/MC);
    RMS_7(:,i) = sqrt(sum(e_x_7(:,i).^2)/MC);
    RMS_8(:,i) = sqrt(sum(e_x_8(:,i).^2)/MC);
    RMS_9(:,i) = sqrt(sum(e_x_9(:,i).^2)/MC);
    RMS_10(:,i) = sqrt(sum(e_x_10(:,i).^2)/MC);
    RMS_11(:,i) = sqrt(sum(e_x_11(:,i).^2)/MC);
    RMS_12(:,i) = sqrt(sum(e_x_12(:,i).^2)/MC);
end

% settings
default_font_size = 12.5; % default font size
default_line_width = 2; % default line width

% (a,b,k) = (1,0,0), 
figure('DefaultAxesFontSize',default_font_size);
subplot(2,1,1)
plot(1:simLen,RMS_1,'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
title('UKF: RMSE $(\alpha, \beta, \kappa) = (1, 0, 0)$','Interpreter','Latex')
subplot(2,1,2)
plot(1:simLen,sgmhat_1,'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
title('UKF: root(P) $(\alpha, \beta, \kappa) = (1, 0, 0)$','Interpreter','Latex')

% (a,b,k) = (0.5,2,0),
figure('DefaultAxesFontSize',default_font_size);
subplot(2,1,1)
plot(1:simLen,RMS_2,'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
title('UKF: RMSE $(\alpha, \beta, \kappa) = (0.5, 2, 0)$','Interpreter','Latex')
subplot(2,1,2)
plot(1:simLen,sgmhat_2,'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
title('UKF: root(P) $(\alpha, \beta, \kappa) = (0.5, 2, 0)$','Interpreter','Latex')

% (a, b, k) = (0, 1, 0)
figure('DefaultAxesFontSize',default_font_size);
subplot(2,1,1)
plot(1:simLen,RMS_3,'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
title('UKF: RMSE $(\alpha, \beta, \kappa) = (1, 0, 1)$','Interpreter','Latex')
subplot(2,1,2)
plot(1:simLen,sgmhat_3,'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
title('UKF: root(P) $(\alpha, \beta, \kappa) = (1, 0, 1)$','Interpreter','Latex')

% (a,b,k) = (0.5,2,1),
figure('DefaultAxesFontSize',default_font_size);
subplot(2,1,1)
plot(1:simLen,RMS_4,'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
title('UKF: RMSE $(\alpha, \beta, \kappa) = (0.5, 2, 1)$','Interpreter','Latex')
subplot(2,1,2)
plot(1:simLen,sgmhat_4,'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
title('UKF: root(P) $(\alpha, \beta, \kappa) = (0.5, 2, 1)$','Interpreter','Latex')

% (a,b,k) = (1,0,2), 
figure('DefaultAxesFontSize',default_font_size);
subplot(2,1,1)
plot(1:simLen,RMS_5,'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
title('UKF: RMSE $(\alpha, \beta, \kappa) = (1, 0, 2)$','Interpreter','Latex')
subplot(2,1,2)
plot(1:simLen,sgmhat_5,'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
title('UKF: root(P) $(\alpha, \beta, \kappa) = (1, 0, 2)$','Interpreter','Latex')

% (a,b,k) = (0.5,2,2),
figure('DefaultAxesFontSize',default_font_size);
subplot(2,1,1)
plot(1:simLen,RMS_6,'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
title('UKF: RMSE $(\alpha, \beta, \kappa) = (0.5, 2, 2)$','Interpreter','Latex')
subplot(2,1,2)
plot(1:simLen,sgmhat_6,'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
title('UKF: root(P) $(\alpha, \beta, \kappa) = (0.5, 2, 2)$','Interpreter','Latex')

% (a, b, k) = (1, 0, 3)
figure('DefaultAxesFontSize',default_font_size);
subplot(2,1,1)
plot(1:simLen,RMS_7,'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
title('UKF: RMSE $(\alpha, \beta, \kappa) = (1, 0, 3)$','Interpreter','Latex')
subplot(2,1,2)
plot(1:simLen,sgmhat_7,'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
title('UKF: root(P) $(\alpha, \beta, \kappa) = (1, 0, 3)$','Interpreter','Latex')

% (a,b,k) = (0.5,2,3),
figure('DefaultAxesFontSize',default_font_size);
subplot(2,1,1)
plot(1:simLen,RMS_4,'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
title('UKF: RMSE $(\alpha, \beta, \kappa) = (0.5, 2, 3)$','Interpreter','Latex')
subplot(2,1,2)
plot(1:simLen,sgmhat_8,'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
title('UKF: root(P) $(\alpha, \beta, \kappa) = (0.5, 2, 3)$','Interpreter','Latex')

% (a,b,k) = (1,0,8), 
figure('DefaultAxesFontSize',default_font_size);
subplot(2,1,1)
plot(1:simLen,RMS_9,'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
title('UKF: RMSE $(\alpha, \beta, \kappa) = (1, 0, 8)$','Interpreter','Latex')
subplot(2,1,2)
plot(1:simLen,sgmhat_9,'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
title('UKF: root(P) $(\alpha, \beta, \kappa) = (1, 0, 8)$','Interpreter','Latex')

% (a,b,k) = (0.5,2,8),
figure('DefaultAxesFontSize',default_font_size);
subplot(2,1,1)
plot(1:simLen,RMS_10,'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
title('UKF: RMSE $(\alpha, \beta, \kappa) = (0.5, 2, 8)$','Interpreter','Latex')
subplot(2,1,2)
plot(1:simLen,sgmhat_10,'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
title('UKF: root(P) $(\alpha, \beta, \kappa) = (0.5, 2, 8)$','Interpreter','Latex')

% (a, b, k) = (1, 0, 15)
figure('DefaultAxesFontSize',default_font_size);
subplot(2,1,1)
plot(1:simLen,RMS_11,'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
title('UKF: RMSE $(\alpha, \beta, \kappa) = (1, 0, 15)$','Interpreter','Latex')
subplot(2,1,2)
plot(1:simLen,sgmhat_11,'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
title('UKF: root(P) $(\alpha, \beta, \kappa) = (1, 0, 15)$','Interpreter','Latex')

% (a,b,k) = (0.5,2,15),
figure('DefaultAxesFontSize',default_font_size);
subplot(2,1,1)
plot(1:simLen,RMS_12,'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
title('UKF: RMSE $(\alpha, \beta, \kappa) = (0.5, 2, 15)$','Interpreter','Latex')
subplot(2,1,2)
plot(1:simLen,sgmhat_12,'LineWidth',default_line_width)
grid on
xlabel('time [sec]')
title('UKF: root(P) $(\alpha, \beta, \kappa) = (0.5, 2, 15)$','Interpreter','Latex')

%% computation time
figure('DefaultAxesFontSize',default_font_size);
plot(1:MC,time_spent_1,'LineWidth',default_line_width)
hold on
plot(1:MC,time_spent_2,'LineWidth',default_line_width)
xlabel('Monte-Carlo iteration step [-]')
ylabel('time [sec]')
legend({'1','2'},'Location','best')
title('computation time')

figure('DefaultAxesFontSize',default_font_size);
plot(1:MC,time_spent_3,'LineWidth',default_line_width)
hold on
plot(1:MC,time_spent_4,'LineWidth',default_line_width)
xlabel('Monte-Carlo iteration step [-]')
ylabel('time [sec]')
legend({'3','4'},'Location','best')
title('computation time')
%%
plot(1:MC,time_spent_5,'LineWidth',default_line_width)
hold on
plot(1:MC,time_spent_6,'LineWidth',default_line_width)
hold on
plot(1:MC,time_spent_7,'LineWidth',default_line_width)
hold on
plot(1:MC,time_spent_8,'LineWidth',default_line_width)
hold on
plot(1:MC,time_spent_9,'LineWidth',default_line_width)
hold on
plot(1:MC,time_spent_10,'LineWidth',default_line_width)
hold on
plot(1:MC,time_spent_11,'LineWidth',default_line_width)
hold on
plot(1:MC,time_spent_12,'LineWidth',default_line_width)
grid on
xlabel('Monte-Carlo iteration step [-]')
ylabel('time [sec]')
legend({'1','2','3','4','5','6','7','8','9','10','11','12'},'Location','best')
title('computation time')
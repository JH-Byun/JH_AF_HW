load 1_f_PF_SIR_1_1.mat
e_x_1_1 = x_true - xhat;
ESS_1_1 = ESS;

load 1_f_PF_SIR_1_2.mat
e_x_1_2 = x_true - xhat;
ESS_1_2 = ESS;

load 1_f_PF_SIR_1_3.mat
e_x_1_3 = x_true - xhat;
ESS_1_3 = ESS;

load 1_f_PF_SIR_1_4.mat
e_x_1_4 = x_true - xhat;
ESS_1_4 = ESS;

load 1_f_PF_SIR_2_1.mat
e_x_2_1 = x_true - xhat;
ESS_2_1 = ESS;

load 1_f_PF_SIR_2_2.mat
e_x_2_2 = x_true - xhat;
ESS_2_2 = ESS;

load 1_f_PF_SIR_2_3.mat
e_x_2_3 = x_true - xhat;
ESS_2_3 = ESS;

load 1_f_PF_SIR_2_4.mat
e_x_2_4 = x_true - xhat;
ESS_2_4 = ESS;

load 1_f_PF_SIR_3_1.mat
e_x_3_1 = x_true - xhat;
ESS_3_1 = ESS;

load 1_f_PF_SIR_3_2.mat
e_x_3_2 = x_true - xhat;
ESS_3_2 = ESS;

load 1_f_PF_SIR_3_3.mat
e_x_3_3 = x_true - xhat;
ESS_3_3 = ESS;

load 1_f_PF_SIR_3_4.mat
e_x_3_4 = x_true - xhat;
ESS_3_4 = ESS;

load 1_f_PF_SIR_4_1.mat
e_x_4_1 = x_true - xhat;
ESS_4_1 = ESS;

load 1_f_PF_SIR_4_2.mat
e_x_4_2 = x_true - xhat;
ESS_4_2 = ESS;

load 1_f_PF_SIR_4_3.mat
e_x_4_3 = x_true - xhat;
ESS_4_3 = ESS;

load 1_f_PF_SIR_4_4.mat
e_x_4_4 = x_true - xhat;
ESS_4_4 = ESS;

MC = size(e_x_1_1,1)/2;
simLen = size(e_x_1_1,2);

for i = 1:simLen
    RMS_SIR_1_1(:,i) = [sqrt(sum(e_x_1_1(1:2:end,i).^2)/MC);sqrt(sum(e_x_1_1(2:2:end,i).^2)/MC)];
    RMS_SIR_1_2(:,i) = [sqrt(sum(e_x_1_2(1:2:end,i).^2)/MC);sqrt(sum(e_x_1_2(2:2:end,i).^2)/MC)];
    RMS_SIR_1_3(:,i) = [sqrt(sum(e_x_1_3(1:2:end,i).^2)/MC);sqrt(sum(e_x_1_3(2:2:end,i).^2)/MC)];
    RMS_SIR_1_4(:,i) = [sqrt(sum(e_x_1_4(1:2:end,i).^2)/MC);sqrt(sum(e_x_1_4(2:2:end,i).^2)/MC)];
    RMS_SIR_2_1(:,i) = [sqrt(sum(e_x_2_1(1:2:end,i).^2)/MC);sqrt(sum(e_x_2_1(2:2:end,i).^2)/MC)];
    RMS_SIR_2_2(:,i) = [sqrt(sum(e_x_2_2(1:2:end,i).^2)/MC);sqrt(sum(e_x_2_2(2:2:end,i).^2)/MC)];
    RMS_SIR_2_3(:,i) = [sqrt(sum(e_x_2_3(1:2:end,i).^2)/MC);sqrt(sum(e_x_2_3(2:2:end,i).^2)/MC)];
    RMS_SIR_2_4(:,i) = [sqrt(sum(e_x_2_4(1:2:end,i).^2)/MC);sqrt(sum(e_x_2_4(2:2:end,i).^2)/MC)];
    RMS_SIR_3_1(:,i) = [sqrt(sum(e_x_3_1(1:2:end,i).^2)/MC);sqrt(sum(e_x_3_1(2:2:end,i).^2)/MC)];
    RMS_SIR_3_2(:,i) = [sqrt(sum(e_x_3_2(1:2:end,i).^2)/MC);sqrt(sum(e_x_3_2(2:2:end,i).^2)/MC)];
    RMS_SIR_3_3(:,i) = [sqrt(sum(e_x_3_3(1:2:end,i).^2)/MC);sqrt(sum(e_x_3_3(2:2:end,i).^2)/MC)];
    RMS_SIR_3_4(:,i) = [sqrt(sum(e_x_3_4(1:2:end,i).^2)/MC);sqrt(sum(e_x_3_4(2:2:end,i).^2)/MC)];
    RMS_SIR_4_1(:,i) = [sqrt(sum(e_x_4_1(1:2:end,i).^2)/MC);sqrt(sum(e_x_4_1(2:2:end,i).^2)/MC)];
    RMS_SIR_4_2(:,i) = [sqrt(sum(e_x_4_2(1:2:end,i).^2)/MC);sqrt(sum(e_x_4_2(2:2:end,i).^2)/MC)];
    RMS_SIR_4_3(:,i) = [sqrt(sum(e_x_4_3(1:2:end,i).^2)/MC);sqrt(sum(e_x_4_3(2:2:end,i).^2)/MC)];
    RMS_SIR_4_4(:,i) = [sqrt(sum(e_x_4_4(1:2:end,i).^2)/MC);sqrt(sum(e_x_4_4(2:2:end,i).^2)/MC)];
end

% settings
default_font_size = 12.5; % default font size
default_line_width = 2; % default line width

% case 1
    % case 1-1
    figure('DefaultAxesFontSize',default_font_size);
    subplot(2,1,1)
    plot(1:simLen,RMS_SIR_1_1(1,:),'LineWidth',default_line_width)
    hold on
    plot(1:simLen,RMS_SIR_1_1(2,:),'LineWidth',default_line_width)
    grid on
    xlabel('time [sec]')
    legend({'RMSE of position','RMSE of velocity'},'Location', 'best')
    subplot(2,1,2)
    plot(1:simLen,ESS_1_1,'LineWidth',default_line_width)
    grid on
    xlabel('time [sec]')
    ylabel('ESS [-]')
    sgtitle('SIR result (q = 10, r = 5, ESS threshold: 100)')
    % case 1-2
    figure('DefaultAxesFontSize',default_font_size);
    subplot(2,1,1)
    plot(1:simLen,RMS_SIR_1_2(1,:),'LineWidth',default_line_width)
    hold on
    plot(1:simLen,RMS_SIR_1_2(2,:),'LineWidth',default_line_width)
    grid on
    xlabel('time [sec]')
    legend({'RMSE of position','RMSE of velocity'},'Location', 'best')
    subplot(2,1,2)
    plot(1:simLen,ESS_1_2,'LineWidth',default_line_width)
    grid on
    xlabel('time [sec]')
    ylabel('ESS [-]')
    sgtitle('SIR result (q = 10, r = 5, ESS threshold: 300)')
    % case 1-3
    figure('DefaultAxesFontSize',default_font_size);
    subplot(2,1,1)
    plot(1:simLen,RMS_SIR_1_3(1,:),'LineWidth',default_line_width)
    hold on
    plot(1:simLen,RMS_SIR_1_3(2,:),'LineWidth',default_line_width)
    grid on
    xlabel('time [sec]')
    legend({'RMSE of position','RMSE of velocity'},'Location', 'best')
    subplot(2,1,2)
    plot(1:simLen,ESS_1_3,'LineWidth',default_line_width)
    grid on
    xlabel('time [sec]')
    ylabel('ESS [-]')
    sgtitle('SIR result (q = 10, r = 5, ESS threshold: 450)')
    % case 1-4
    figure('DefaultAxesFontSize',default_font_size);
    subplot(2,1,1)
    plot(1:simLen,RMS_SIR_1_4(1,:),'LineWidth',default_line_width)
    hold on
    plot(1:simLen,RMS_SIR_1_4(2,:),'LineWidth',default_line_width)
    grid on
    xlabel('time [sec]')
    legend({'RMSE of position','RMSE of velocity'},'Location', 'best')
    subplot(2,1,2)
    plot(1:simLen,ESS_1_4,'LineWidth',default_line_width)
    grid on
    xlabel('time [sec]')
    ylabel('ESS [-]')
    sgtitle('SIR result (q = 10, r = 5, ESS threshold: 600)')

% case 2
    % case 2-1
    figure('DefaultAxesFontSize',default_font_size);
    subplot(2,1,1)
    plot(1:simLen,RMS_SIR_2_1(1,:),'LineWidth',default_line_width)
    hold on
    plot(1:simLen,RMS_SIR_2_1(2,:),'LineWidth',default_line_width)
    grid on
    xlabel('time [sec]')
    legend({'RMSE of position','RMSE of velocity'},'Location', 'best')
    subplot(2,1,2)
    plot(1:simLen,ESS_2_1,'LineWidth',default_line_width)
    grid on
    xlabel('time [sec]')
    ylabel('ESS [-]')
    sgtitle('SIR result (q = 1, r = 5, ESS threshold: 100)')
    % case 2-2
    figure('DefaultAxesFontSize',default_font_size);
    subplot(2,1,1)
    plot(1:simLen,RMS_SIR_2_2(1,:),'LineWidth',default_line_width)
    hold on
    plot(1:simLen,RMS_SIR_2_2(2,:),'LineWidth',default_line_width)
    grid on
    xlabel('time [sec]')
    legend({'RMSE of position','RMSE of velocity'},'Location', 'best')
    subplot(2,1,2)
    plot(1:simLen,ESS_2_2,'LineWidth',default_line_width)
    grid on
    xlabel('time [sec]')
    ylabel('ESS [-]')
    sgtitle('SIR result (q = 1, r = 5, ESS threshold: 300)')
    % case 2-3
    figure('DefaultAxesFontSize',default_font_size);
    subplot(2,1,1)
    plot(1:simLen,RMS_SIR_2_3(1,:),'LineWidth',default_line_width)
    hold on
    plot(1:simLen,RMS_SIR_2_3(2,:),'LineWidth',default_line_width)
    grid on
    xlabel('time [sec]')
    legend({'RMSE of position','RMSE of velocity'},'Location', 'best')
    subplot(2,1,2)
    plot(1:simLen,ESS_2_3,'LineWidth',default_line_width)
    grid on
    xlabel('time [sec]')
    ylabel('ESS [-]')
    sgtitle('SIR result (q = 1, r = 5, ESS threshold: 450)')
    % case 2-4
    figure('DefaultAxesFontSize',default_font_size);
    subplot(2,1,1)
    plot(1:simLen,RMS_SIR_2_4(1,:),'LineWidth',default_line_width)
    hold on
    plot(1:simLen,RMS_SIR_2_4(2,:),'LineWidth',default_line_width)
    grid on
    xlabel('time [sec]')
    legend({'RMSE of position','RMSE of velocity'},'Location', 'best')
    subplot(2,1,2)
    plot(1:simLen,ESS_2_4,'LineWidth',default_line_width)
    grid on
    xlabel('time [sec]')
    ylabel('ESS [-]')
    sgtitle('SIR result (q = 1, r = 5, ESS threshold: 600)')    
    
% case 3
    % case 3-1
    figure('DefaultAxesFontSize',default_font_size);
    subplot(2,1,1)
    plot(1:simLen,RMS_SIR_3_1(1,:),'LineWidth',default_line_width)
    hold on
    plot(1:simLen,RMS_SIR_3_1(2,:),'LineWidth',default_line_width)
    grid on
    xlabel('time [sec]')
    legend({'RMSE of position','RMSE of velocity'},'Location', 'best')
    subplot(2,1,2)
    plot(1:simLen,ESS_3_1,'LineWidth',default_line_width)
    grid on
    xlabel('time [sec]')
    ylabel('ESS [-]')
    sgtitle('SIR result (q = 0.1, r = 5, ESS threshold: 100)')
    % case 3-2
    figure('DefaultAxesFontSize',default_font_size);
    subplot(2,1,1)
    plot(1:simLen,RMS_SIR_3_2(1,:),'LineWidth',default_line_width)
    hold on
    plot(1:simLen,RMS_SIR_3_2(2,:),'LineWidth',default_line_width)
    grid on
    xlabel('time [sec]')
    legend({'RMSE of position','RMSE of velocity'},'Location', 'best')
    subplot(2,1,2)
    plot(1:simLen,ESS_3_2,'LineWidth',default_line_width)
    grid on
    xlabel('time [sec]')
    ylabel('ESS [-]')
    sgtitle('SIR result (q = 0.1, r = 5, ESS threshold: 300)')
    % case 3-3
    figure('DefaultAxesFontSize',default_font_size);
    subplot(2,1,1)
    plot(1:simLen,RMS_SIR_3_3(1,:),'LineWidth',default_line_width)
    hold on
    plot(1:simLen,RMS_SIR_3_3(2,:),'LineWidth',default_line_width)
    grid on
    xlabel('time [sec]')
    legend({'RMSE of position','RMSE of velocity'},'Location', 'best')
    subplot(2,1,2)
    plot(1:simLen,ESS_3_3,'LineWidth',default_line_width)
    grid on
    xlabel('time [sec]')
    ylabel('ESS [-]')
    sgtitle('SIR result (q = 0.1, r = 5, ESS threshold: 450)')
    % case 3-4
    figure('DefaultAxesFontSize',default_font_size);
    subplot(2,1,1)
    plot(1:simLen,RMS_SIR_3_4(1,:),'LineWidth',default_line_width)
    hold on
    plot(1:simLen,RMS_SIR_3_4(2,:),'LineWidth',default_line_width)
    grid on
    xlabel('time [sec]')
    legend({'RMSE of position','RMSE of velocity'},'Location', 'best')
    subplot(2,1,2)
    plot(1:simLen,ESS_3_4,'LineWidth',default_line_width)
    grid on
    xlabel('time [sec]')
    ylabel('ESS [-]')
    sgtitle('SIR result (q = 0.1, r = 5, ESS threshold: 600)')    
    
% case 4
    % case 4-1
    figure('DefaultAxesFontSize',default_font_size);
    subplot(2,1,1)
    plot(1:simLen,RMS_SIR_4_1(1,:),'LineWidth',default_line_width)
    hold on
    plot(1:simLen,RMS_SIR_4_1(2,:),'LineWidth',default_line_width)
    grid on
    xlabel('time [sec]')
    legend({'RMSE of position','RMSE of velocity'},'Location', 'best')
    subplot(2,1,2)
    plot(1:simLen,ESS_4_1,'LineWidth',default_line_width)
    grid on
    xlabel('time [sec]')
    ylabel('ESS [-]')
    sgtitle('SIR result (q = 0.01, r = 5, ESS threshold: 100)')
    % case 4-2
    figure('DefaultAxesFontSize',default_font_size);
    subplot(2,1,1)
    plot(1:simLen,RMS_SIR_4_2(1,:),'LineWidth',default_line_width)
    hold on
    plot(1:simLen,RMS_SIR_4_2(2,:),'LineWidth',default_line_width)
    grid on
    xlabel('time [sec]')
    legend({'RMSE of position','RMSE of velocity'},'Location', 'best')
    subplot(2,1,2)
    plot(1:simLen,ESS_4_2,'LineWidth',default_line_width)
    grid on
    xlabel('time [sec]')
    ylabel('ESS [-]')
    sgtitle('SIR result (q = 0.01, r = 5, ESS threshold: 300)')
    % case 4-3
    figure('DefaultAxesFontSize',default_font_size);
    subplot(2,1,1)
    plot(1:simLen,RMS_SIR_4_3(1,:),'LineWidth',default_line_width)
    hold on
    plot(1:simLen,RMS_SIR_4_3(2,:),'LineWidth',default_line_width)
    grid on
    xlabel('time [sec]')
    legend({'RMSE of position','RMSE of velocity'},'Location', 'best')
    subplot(2,1,2)
    plot(1:simLen,ESS_4_3,'LineWidth',default_line_width)
    grid on
    xlabel('time [sec]')
    ylabel('ESS [-]')
    sgtitle('SIR result (q = 0.01, r = 5, ESS threshold: 450)')
    % case 4-4
    figure('DefaultAxesFontSize',default_font_size);
    subplot(2,1,1)
    plot(1:simLen,RMS_SIR_4_4(1,:),'LineWidth',default_line_width)
    hold on
    plot(1:simLen,RMS_SIR_4_4(2,:),'LineWidth',default_line_width)
    grid on
    xlabel('time [sec]')
    legend({'RMSE of position','RMSE of velocity'},'Location', 'best')
    subplot(2,1,2)
    plot(1:simLen,ESS_4_4,'LineWidth',default_line_width)
    grid on
    xlabel('time [sec]')
    ylabel('ESS [-]')
    sgtitle('SIR result (q = 0.01, r = 5, ESS threshold: 600)')    
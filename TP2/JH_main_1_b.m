clear all; close all; clc

addpath('./resampling methods from Tiancheng Li');
addpath('./ekfukf');

%% Parameters
MC = 100;
simLen = 50;
x0 = [0;0];
sqrtP0 = [5 0;...
    0 1];
T = 1;
F = [1 T;...
    0 1];
Gamma = [T^2/2;T];
H = [1 0];
sqrtQ = 1;
sqrtR = 5;
measMC = 1;

%% simulation (ensemble generation)
x_true = zeros(MC*2, simLen);
z_true = zeros(MC, simLen);
z = zeros(MC, simLen);
for mcid = 1:MC
    rng(mcid);
    x_true(mcid*2-1:mcid*2,1) = x0 + sqrtP0*randn(2,1);
    for i=2:simLen
        x_true(mcid*2-1:mcid*2,i) = F*x_true(mcid*2-1:mcid*2,i-1) + Gamma*sqrtQ*randn(1,1);
    end
    z_true(mcid,:) = x_true(mcid*2-1,:);
    z(measMC*(mcid-1)+1:measMC*mcid,:) = repmat(z_true(mcid,:),measMC,1) + sqrtR*randn(measMC,simLen);
end

%% filtering
xtrue_array = reshape(repmat(reshape(x_true,2,[]),measMC,1),MC*measMC*2,[]);

% UKF
xhat_UKF = zeros(MC*2, simLen);
sgmhat_UKF = zeros(MC*2, simLen);
time_spent_UKF = zeros(1, MC);

for mcid = 1:MC
    tic;
    xhat_UKF(2*mcid - 1:2*mcid,1) = x0;
    Phat = sqrtP0^2;
    for i = 2:simLen
        [xhat_UKF(2*mcid - 1:2*mcid,i), Phat] = ukf_predict1(xhat_UKF(2*mcid - 1:2*mcid,i-1), Phat, F, Gamma*sqrtQ^2*Gamma.',[],0.5,2,1,[]);
        [xhat_UKF(2*mcid - 1:2*mcid,i), Phat] = ukf_update1(xhat_UKF(2*mcid - 1:2*mcid,i), Phat, z(mcid, i), [1 0], sqrtR^2, [], 0.5, 2, 1, []);
        sgmhat_UKF(2*mcid - 1:2*mcid,i) = [sqrt(Phat(1,1));sqrt(Phat(2,2))];
    end
    time_spent_UKF(1,mcid) = toc;
end
fprintf('UKF done ... \n');

% SIS 
xhat_SIS = zeros(MC*2, simLen);
sgmhat_SIS = zeros(MC*2, simLen);
time_spent_SIS = zeros(1, MC);
N = 1000; % particle numbers

for mcid = 1:MC
    tic;
    wgt = 1/N*ones(1,N);
    xpart_prior = gauss_rnd(x0,sqrtP0^2,N);
    xhat_SIS(2*mcid-1:2*mcid,1) = sum(xpart_prior.*repmat(wgt,2,1),2);
    sgmhat_SIS(2*mcid-1:2*mcid,1) = [sqrtP0(1,1);sqrtP0(2,2)];
    for i = 2:simLen
        xpart_prior = F*xpart_prior + Gamma*gauss_rnd(0,sqrtQ^2,N);
        ypart = H * xpart_prior;
        for j = 1:N
            p_yx(1,j) = gauss_pdf(z(mcid, i),ypart(1,j),sqrtR^2).';
        end
        wgt = p_yx.*wgt;
        wgt = wgt/N;
        wgt = wgt/sum(wgt);
        
        xhat_SIS(2*mcid-1:2*mcid,i) = sum(xpart_prior.*repmat(wgt,2,1),2);
        P = (xpart_prior - repmat(xhat_SIS(2*mcid-1:2*mcid,i),1,N))*diag(wgt)*...
            (xpart_prior - repmat(xhat_SIS(2*mcid-1:2*mcid,i),1,N)).';
        sgmhat_SIS(2*mcid-1,i) = sqrt(P(1,1));
        sgmhat_SIS(2*mcid,i) = sqrt(P(2,2));
    end
    time_spent_SIS(1,mcid) = toc;
end
fprintf('PF - SIS done ... \n');

% SIR 
xhat_SIR = zeros(MC*2, simLen);
sgmhat_SIR = zeros(MC*2, simLen);
time_spent_SIR = zeros(1, MC);
N = 1000; % particle numbers

for mcid = 1:MC
    tic;
    wgt = 1/N*ones(1,N);
    xpart_prior = gauss_rnd(x0,sqrtP0^2,N);
    xhat_SIR(2*mcid-1:2*mcid,1) = sum(xpart_prior.*repmat(wgt,2,1),2);
    sgmhat_SIR(2*mcid-1:2*mcid,1) = [sqrtP0(1,1);sqrtP0(2,2)];
    for i = 2:simLen
        resam_id = resampleSystematic(wgt);
        xpart_prior = xpart_prior(:,resam_id);
        xpart_prior = F*xpart_prior + Gamma*gauss_rnd(0,sqrtQ^2,N);
        
        ypart = H * xpart_prior;
        for j = 1:N
            p_yx(1,j) = gauss_pdf(z(mcid, i),ypart(1,j),sqrtR^2).';
        end
        p_yx = p_yx/sum(p_yx);
        wgt = p_yx/N;
        wgt = wgt/sum(wgt);
        
        xhat_SIR(2*mcid-1:2*mcid,i) = sum(xpart_prior.*repmat(wgt,2,1),2);
        P = (xpart_prior - repmat(xhat_SIR(2*mcid-1:2*mcid,i),1,N))*diag(wgt)*...
            (xpart_prior - repmat(xhat_SIR(2*mcid-1:2*mcid,i),1,N)).';
        sgmhat_SIR(2*mcid-1,i) = sqrt(P(1,1));
        sgmhat_SIR(2*mcid,i) = sqrt(P(2,2));
    end
    time_spent_SIR(1,mcid) = toc;
end
fprintf('PF - SIR done ... \n');

%% results
% settings
load KF.mat
default_font_size = 12.5; % default font size
default_line_width = 2; % default line width
     
% RMS error comparison
    % RMS error
    xerr_UKF = xtrue_array - xhat_UKF;
    xerr_SIS = xtrue_array - xhat_SIS;
    xerr_SIR = xtrue_array - xhat_SIR;     
    for i = 1:simLen
        RMS_UKF(:,i) = [sqrt(sum(xerr_UKF(1:2:end,i).^2)/MC);sqrt(sum(xerr_UKF(2:2:end,i).^2)/MC)];
        RMS_SIS(:,i) = [sqrt(sum(xerr_SIS(1:2:end,i).^2)/MC);sqrt(sum(xerr_SIS(2:2:end,i).^2)/MC)];
        RMS_SIR(:,i) = [sqrt(sum(xerr_SIR(1:2:end,i).^2)/MC);sqrt(sum(xerr_SIR(2:2:end,i).^2)/MC)];
    end
    % RMSE comparison plot
    figure('DefaultAxesFontSize',default_font_size);
    subplot(2,1,1)
    plot(1:simLen,RMS_KF(1,:),'LineWidth',default_line_width)
    hold on
    plot(1:simLen,RMS_UKF(1,:),'LineWidth',default_line_width)
    hold on
    plot(1:simLen,RMS_SIS(1,:),'LineWidth',default_line_width)
    hold on
    plot(1:simLen,RMS_SIR(1,:),'LineWidth',default_line_width)
    grid on
    ylabel('position [m]')
    title('Position RMS error')
    legend({'KF','UKF','PF - SIS','PF - SIR'},'Location','best')
    subplot(2,1,2)
    plot(1:simLen,RMS_KF(2,:),'LineWidth',default_line_width)
    hold on
    plot(1:simLen,RMS_UKF(2,:),'LineWidth',default_line_width)
    hold on
    plot(1:simLen,RMS_SIS(2,:),'LineWidth',default_line_width)
    hold on
    plot(1:simLen,RMS_SIR(2,:),'LineWidth',default_line_width)
    grid on
    xlabel('time [sec]')
    ylabel('velocity [m/s]')
    title('Velocity RMS error')
    
% standard deviation error comparison
    % KF
    figure('DefaultAxesFontSize',default_font_size);
    subplot(2,1,1)
    plot(1:simLen,sgmhat_KF(1:2:end, :),'LineWidth',default_line_width)
    grid on
    ylabel('$\sigma_x$ [m]','Interpreter','Latex')
    title('Position standard deviation (KF)')
    subplot(2,1,2)
    plot(1:simLen,sgmhat_KF(2:2:end, :),'LineWidth',default_line_width)
    grid on
    xlabel('time [sec]')
    ylabel('$\sigma_v$ [m/s]','Interpreter','Latex')
    ylim([0 5])
    title('Velocity standard deviation (KF)')
    
    % UKF
    figure('DefaultAxesFontSize',default_font_size);
    subplot(2,1,1)
    plot(1:simLen,sgmhat_UKF(1:2:end, :),'LineWidth',default_line_width)
    grid on
    ylabel('$\sigma_x$ [m]','Interpreter','Latex')
    title('Position standard deviation (UKF)')
    subplot(2,1,2)
    plot(1:simLen,sgmhat_UKF(2:2:end, :),'LineWidth',default_line_width)
    grid on
    xlabel('time [sec]')
    ylabel('$\sigma_v$ [m/s]','Interpreter','Latex')
    ylim([0 5])
    title('Velocity standard deviation (UKF)')
    
    % PF - SIS
    figure('DefaultAxesFontSize',default_font_size);
    subplot(2,1,1)
    plot(1:simLen,sgmhat_SIS(1:2:end, :),'LineWidth',default_line_width)
    grid on
    ylabel('$\sigma_x$ [m]','Interpreter','Latex')
    title('Position standard deviation (PF - SIS)')
    subplot(2,1,2)
    plot(1:simLen,sgmhat_SIS(2:2:end, :),'LineWidth',default_line_width)
    grid on
    xlabel('time [sec]')
    ylabel('$\sigma_v$ [m/s]','Interpreter','Latex')
    ylim([0 5])
    title('Velocity standard deviation (PF - SIS)')
    
    % PF - SIR
    figure('DefaultAxesFontSize',default_font_size);
    subplot(2,1,1)
    plot(1:simLen,sgmhat_SIR(1:2:end, :),'LineWidth',default_line_width)
    grid on
    ylabel('$\sigma_x$ [m]','Interpreter','Latex')
    title('Position standard deviation (PF - SIR)')
    subplot(2,1,2)
    plot(1:simLen,sgmhat_SIR(2:2:end, :),'LineWidth',default_line_width)
    grid on
    xlabel('time [sec]')
    ylabel('$\sigma_v$ [m/s]','Interpreter','Latex')
    ylim([0 5])
    title('Velocity standard deviation (PF - SIR)')
    
% computation time comparison
    figure('DefaultAxesFontSize',default_font_size);
    plot(1:MC,time_spent_KF(1, :),'LineWidth',default_line_width)
    hold on
    plot(1:MC,time_spent_UKF(1, :),'LineWidth',default_line_width)
    hold on
    plot(1:MC,time_spent_SIS(1, :),'LineWidth',default_line_width)
    hold on
    plot(1:MC,time_spent_SIR(1, :),'LineWidth',default_line_width)
    grid on
    xlabel('MC iteration number [-]')
    ylabel('computation time [sec]')
    legend({'KF','UKF','PF - SIS','PF - SIR'},'Location','best')
    title('Computation time comparison')

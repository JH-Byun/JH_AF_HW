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
        [xhat_UKF(2*mcid - 1:2*mcid,i), Phat] = ukf_predict1(xhat_UKF(2*mcid - 1:2*mcid,i-1), Phat, F, Gamma*sqrtQ^2*Gamma.',[],[],[],1,[]);
        [xhat_UKF(2*mcid - 1:2*mcid,i), Phat] = ukf_update1(xhat_UKF(2*mcid - 1:2*mcid,i), Phat, z(mcid, i), [1 0], sqrtR^2, [], [], [], 1, []);
        sgmhat_UKF(2*mcid - 1:2*mcid,i) = [sqrt(Phat(1,1));sqrt(Phat(2,2))];
    end
    time_spent_UKF(1,mcid) = toc;
end

% SIS 
xhat_SIS = zeros(MC*2, simLen);
sgmhat_SIS = zeros(MC*2, simLen);
time_spent_SIS = zeros(1, MC);

for mcid = 1:MC
    tic;
    for i = 2:simLen
        xpart_prior(:,i) = F*xpart_prior(:,i) + Gamma*sqrtQ*randn(1,1);
        ypart = H*xpart_prior(:,i) + sqrtR*randn(1,1);
        vhat = z(mcid, i) - ypart;
        q(i) = exp(-vhat^2/(2*sqrtR^2))/(sqrtR*sqrt(2*pi));
        wgt(i) = % compute weights 
    end
    time_spent_SIS(1,mcid) = toc;
end

%% results
% settings
default_font_size = 12.5; % default font size
default_line_width = 2; % default line width

% KF
xerr = xtrue_array - xhat;
    % X(1): position
    figure('DefaultAxesFontSize',default_font_size);
    for mcid =  1:MC
        hold on
        plot(1:simLen,xerr(2*mcid - 1, :),'LineWidth',default_line_width)
    end
    grid on
    xlabel('time [sec]')
    ylabel('position [m]')
    title('Position error')
    % X(2): velocity
    figure('DefaultAxesFontSize',default_font_size);
    for mcid =  1:MC
        hold on
        plot(1:simLen,xerr(2*mcid, :),'LineWidth',default_line_width)
    end
    grid on
    xlabel('time [sec]')
    ylabel('velocity [m/s]')
    title('Velocity error')
    
    % RMS error
    RMS = zeros(2,simLen);
    for i = 1:simLen
        S_sum = zeros(2,1);
        for mcid = 1:MC
            S_sum(1) = S_sum(1) + xerr(2*mcid - 1,i)^2;
            S_sum(2) = S_sum(2) + xerr(2*mcid,i)^2;
        end
        RMS(1,i) = sqrt(S_sum(1)/MC);
        RMS(2,i) = sqrt(S_sum(2)/MC);
    end
        % RMS error plot
        figure('DefaultAxesFontSize',default_font_size);
        subplot(2,1,1)
        plot(1:simLen,RMS(1,:),'LineWidth',default_line_width)
        grid on
        xlabel('time [sec]')
        ylabel('position [m]')
        ylim([0 5])
        title('Position RMS error')
        subplot(2,1,2)
        plot(1:simLen,RMS(2,:),'LineWidth',default_line_width)
        grid on
        xlabel('time [sec]')
        ylabel('velocity [m/2]')
        ylim([0 5])
        title('Velocity RMS error')
        
    % sqrtP diagonal term plot
        % sgm_x: position standard deviation
        figure('DefaultAxesFontSize',default_font_size);
        subplot(2,1,1)
        for mcid =  1:MC
            hold on
            plot(1:simLen,sgmhat(2*mcid - 1, :),'LineWidth',default_line_width)
        end
        grid on
        xlabel('time [sec]')
        ylabel('$\sigma_x$ [m]','Interpreter','Latex')
        ylim([0 5])
        title('Position standard deviation')
        % sgm_v: velocity standard deviation
        subplot(2,1,2)
        for mcid =  1:MC
            hold on
            plot(1:simLen,sgmhat(2*mcid, :),'LineWidth',default_line_width)
        end
        grid on
        xlabel('time [sec]')
        ylabel('$\sigma_v$ [m/s]','Interpreter','Latex')
        ylim([0 5])
        title('Velocity standard deviation')
       
 % UKF
xerr_UKF = xtrue_array - xhat_UKF;
    % RMS error
    RMS_UKF = zeros(2,simLen);
    for i = 1:simLen
        S_sum = zeros(2,1);
        for mcid = 1:MC
            S_sum(1) = S_sum(1) + xerr_UKF(2*mcid - 1,i)^2;
            S_sum(2) = S_sum(2) + xerr_UKF(2*mcid,i)^2;
        end
        RMS_UKF(1,i) = sqrt(S_sum(1)/MC);
        RMS_UKF(2,i) = sqrt(S_sum(2)/MC);
    end
        % RMS error plot
        figure('DefaultAxesFontSize',default_font_size);
        subplot(2,1,1)
        plot(1:simLen,RMS_UKF(1,:),'LineWidth',default_line_width)
        grid on
        xlabel('time [sec]')
        ylabel('position [m]')
        ylim([0 5])
        title('Position RMS error (UKF)')
        subplot(2,1,2)
        plot(1:simLen,RMS_UKF(2,:),'LineWidth',default_line_width)
        grid on
        xlabel('time [sec]')
        ylabel('velocity [m/2]')
        ylim([0 5])
        title('Velocity RMS error (UKF)')
        
    % sqrtP diagonal term plot
        % sgm_x: position standard deviation
        figure('DefaultAxesFontSize',default_font_size);
        subplot(2,1,1)
        for mcid =  1:MC
            hold on
            plot(1:simLen,sgmhat_UKF(2*mcid - 1, :),'LineWidth',default_line_width)
        end
        grid on
        xlabel('time [sec]')
        ylabel('$\sigma_x$ [m]','Interpreter','Latex')
        ylim([0 5])
        title('Position standard deviation (UKF)')
        % sgm_v: velocity standard deviation
        subplot(2,1,2)
        for mcid =  1:MC
            hold on
            plot(1:simLen,sgmhat_UKF(2*mcid, :),'LineWidth',default_line_width)
        end
        grid on
        xlabel('time [sec]')
        ylabel('$\sigma_v$ [m/s]','Interpreter','Latex')
        ylim([0 5])
        title('Velocity standard deviation (UKF)')
  
        
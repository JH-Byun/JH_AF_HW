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

% KF
xhat_KF = zeros(MC*2, simLen);
sgmhat_KF = zeros(MC*2, simLen);
time_spent_KF = zeros(1, MC);

for mcid = 1:MC
    tic;
    [xhat_KF(2*mcid - 1:2*mcid,:), Phat] = kf_loop(x0, sqrtP0^2, [1 0], sqrtR^2, z(mcid,:), F, Gamma*sqrtQ^2*Gamma.');
    for i = 1:simLen
        sgmhat_KF(2*mcid - 1:2*mcid,i) = [sqrt(Phat(1,1,i));sqrt(Phat(2,2,i))];
    end
    time_spent_KF(1,mcid) = toc;
end

%% results
% settings
default_font_size = 12.5; % default font size
default_line_width = 2; % default line width

% KF
xerr = xtrue_array - xhat_KF;
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
    for i = 1:simLen
        RMS_KF(:,i) = [sqrt(sum(xerr(1:2:end,i).^2)/MC);sqrt(sum(xerr(2:2:end,i).^2)/MC)];
    end
    % RMS error plot
    figure('DefaultAxesFontSize',default_font_size);
    subplot(2,1,1)
    plot(1:simLen,RMS_KF(1,:),'LineWidth',default_line_width)
    grid on
    ylabel('position [m]')
    ylim([0 5])
    title('Position RMS error')
    subplot(2,1,2)
    plot(1:simLen,RMS_KF(2,:),'LineWidth',default_line_width)
    grid on
    xlabel('time [sec]')
    ylabel('velocity [m/s]')
    ylim([0 5])
    title('Velocity RMS error')
    
    figure('DefaultAxesFontSize',default_font_size);
    subplot(2,1,1)
    plot(1:simLen,sgmhat_KF(1:2:end, :),'LineWidth',default_line_width)
    grid on
    xlabel('time [sec]')
    ylabel('$\sigma_x$ [m]','Interpreter','Latex')
    ylim([0 5])
    title('Position standard deviation')
    subplot(2,1,2)
    plot(1:simLen,sgmhat_KF(2:2:end, :),'LineWidth',default_line_width)
    grid on
    xlabel('time [sec]')
    ylabel('$\sigma_v$ [m/s]','Interpreter','Latex')
    ylim([0 5])
    title('Velocity standard deviation')
    
%% Data acquisition
save KF.mat RMS_KF sgmhat_KF time_spent_KF
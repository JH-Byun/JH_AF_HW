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
N = 1000; % particle number

%% simulation (ensemble generation) + filtering
% case 1 - sqrtQ = 10, sqrtR = 5
sqrtQ = 10;
sqrtR = 5;
    % simulation
    x_true_1 = zeros(2, simLen);
    z_true_1 = zeros(1, simLen);
    z_1 = zeros(1, simLen);
    x_true_1(:,1) = x0 + sqrtP0*randn(2,1);
    for i=2:simLen
        x_true_1(:,i) = F*x_true_1(:,i-1) + Gamma*sqrtQ*randn(1,1);
    end
    z_true_1(1,:) = x_true_1(1,:);
    z_1(1,:) = z_true_1 + sqrtR*randn(1,simLen);
    
    % filtering (SIS)
    xhat_SIS_1 = zeros(2, simLen);
    sgmhat_SIS_1 = zeros(2, simLen);

    tic;
    wgt = 1/N*ones(1,N);
    xpart_prior = gauss_rnd(x0,sqrtP0^2,N);
    xhat_SIS_1(:,1) = sum(xpart_prior,2);
    sgmhat_SIS_1(:,1) = [sqrtP0(1,1);sqrtP0(2,2)];
    for i = 2:simLen
        xpart_prior = F*xpart_prior + Gamma*gauss_rnd(0,sqrtQ^2,N);
        ypart = H * xpart_prior;
        for j = 1:N
            p_yx(1,j) = gauss_pdf(z_1(1, i),ypart(1,j),sqrtR^2).';
        end
        wgt = p_yx.*wgt;
        wgt = wgt/N;
        wgt = wgt/sum(wgt);

        xhat_SIS_1(:,i) = sum(xpart_prior.*repmat(wgt,2,1),2);
        P = (xpart_prior - repmat(xhat_SIS_1(:,i),1,N))*diag(wgt)*...
            (xpart_prior - repmat(xhat_SIS_1(:,i),1,N)).';
        sgmhat_SIS_1(1,i) = sqrt(P(1,1));
        sgmhat_SIS_1(2,i) = sqrt(P(2,2));
    end
    time_spent_SIS_1 = toc;
    
    % filtering (SIR)
    xhat_SIR_1 = zeros(2, simLen);
    sgmhat_SIR_1 = zeros(2, simLen);

    tic;
    wgt = 1/N*ones(1,N);
    xpart_prior = gauss_rnd(x0,sqrtP0^2,N);
    xhat_SIR_1(:,1) = sum(xpart_prior,2);
    sgmhat_SIR_1(:,1) = [sqrtP0(1,1);sqrtP0(2,2)];
    for i = 2:simLen
        resam_id = resampleSystematic(wgt);
        xpart_prior = xpart_prior(:,resam_id);
        xpart_prior = F*xpart_prior + Gamma*gauss_rnd(0,sqrtQ^2,N);

        ypart = H * xpart_prior;
        for j = 1:N
            p_yx(1,j) = gauss_pdf(z_1(1, i),ypart(1,j),sqrtR^2).';
        end
        wgt = p_yx.*wgt;
        wgt = wgt/N;
        wgt = wgt/sum(wgt);

        xhat_SIR_1(:,i) = sum(xpart_prior.*repmat(wgt,2,1),2);
        P = (xpart_prior - repmat(xhat_SIR_1(:,i),1,N))*diag(wgt)*...
            (xpart_prior - repmat(xhat_SIR_1(:,i),1,N)).';
        sgmhat_SIR_1(1,i) = sqrt(P(1,1));
        sgmhat_SIR_1(2,i) = sqrt(P(2,2));
    end
    time_spent_SIR_1 = toc;
    
% case 2 - sqrtQ = 1, sqrtR = 5
sqrtQ = 1;
sqrtR = 5;
    % simulation
    x_true_2 = zeros(2, simLen);
    z_true_2 = zeros(1, simLen);
    z_2 = zeros(1, simLen);
    x_true_2(:,1) = x0 + sqrtP0*randn(2,1);
    for i=2:simLen
        x_true_2(:,i) = F*x_true_2(:,i-1) + Gamma*sqrtQ*randn(1,1);
    end
    z_true_2(1,:) = x_true_2(1,:);
    z_2(1,:) = z_true_2 + sqrtR*randn(1,simLen);
    
    % filtering (SIS)
    xhat_SIS_2 = zeros(2, simLen);
    sgmhat_SIS_2 = zeros(2, simLen);

    tic;
    wgt = 1/N*ones(1,N);
    xpart_prior = gauss_rnd(x0,sqrtP0^2,N);
    xhat_SIS_2(:,1) = sum(xpart_prior,2);
    sgmhat_SIS_2(:,1) = [sqrtP0(1,1);sqrtP0(2,2)];
    for i = 2:simLen
        xpart_prior = F*xpart_prior + Gamma*gauss_rnd(0,sqrtQ^2,N);
        ypart = H * xpart_prior;
        for j = 1:N
            p_yx(1,j) = gauss_pdf(z_2(1, i),ypart(1,j),sqrtR^2).';
        end
        wgt = p_yx.*wgt;
        wgt = wgt/N;
        wgt = wgt/sum(wgt);

        xhat_SIS_2(:,i) = sum(xpart_prior.*repmat(wgt,2,1),2);
        P = (xpart_prior - repmat(xhat_SIS_2(:,i),1,N))*diag(wgt)*...
            (xpart_prior - repmat(xhat_SIS_2(:,i),1,N)).';
        sgmhat_SIS_2(1,i) = sqrt(P(1,1));
        sgmhat_SIS_2(2,i) = sqrt(P(2,2));
    end
    time_spent_SIS_2 = toc;
    
    % filtering (SIR)
    xhat_SIR_2 = zeros(2, simLen);
    sgmhat_SIR_2 = zeros(2, simLen);

    tic;
    wgt = 1/N*ones(1,N);
    xpart_prior = gauss_rnd(x0,sqrtP0^2,N);
    xhat_SIR_2(:,1) = sum(xpart_prior,2);
    sgmhat_SIR_2(:,1) = [sqrtP0(1,1);sqrtP0(2,2)];
    for i = 2:simLen
        resam_id = resampleSystematic(wgt);
        xpart_prior = xpart_prior(:,resam_id);
        xpart_prior = F*xpart_prior + Gamma*gauss_rnd(0,sqrtQ^2,N);

        ypart = H * xpart_prior;
        for j = 1:N
            p_yx(1,j) = gauss_pdf(z_2(1, i),ypart(1,j),sqrtR^2).';
        end
        wgt = p_yx.*wgt;
        wgt = wgt/N;
        wgt = wgt/sum(wgt);

        xhat_SIR_2(:,i) = sum(xpart_prior.*repmat(wgt,2,1),2);
        P = (xpart_prior - repmat(xhat_SIR_2(:,i),1,N))*diag(wgt)*...
            (xpart_prior - repmat(xhat_SIR_2(:,i),1,N)).';
        sgmhat_SIR_2(1,i) = sqrt(P(1,1));
        sgmhat_SIR_2(2,i) = sqrt(P(2,2));
    end
    time_spent_SIR_2 = toc;
    
% case 3 - sqrtQ = 0.1, sqrtR = 5
sqrtQ = 0.1;
sqrtR = 5;
    % simulation
    x_true_3 = zeros(2, simLen);
    z_true_3 = zeros(1, simLen);
    z_3 = zeros(1, simLen);
    x_true_3(:,1) = x0 + sqrtP0*randn(2,1);
    for i=2:simLen
        x_true_3(:,i) = F*x_true_3(:,i-1) + Gamma*sqrtQ*randn(1,1);
    end
    z_true_3(1,:) = x_true_3(1,:);
    z_3(1,:) = z_true_3 + sqrtR*randn(1,simLen);
    
    % filtering (SIS)
    xhat_SIS_3 = zeros(2, simLen);
    sgmhat_SIS_3 = zeros(2, simLen);

    tic;
    wgt = 1/N*ones(1,N);
    xpart_prior = gauss_rnd(x0,sqrtP0^2,N);
    xhat_SIS_3(:,1) = sum(xpart_prior,2);
    sgmhat_SIS_3(:,1) = [sqrtP0(1,1);sqrtP0(2,2)];
    for i = 2:simLen
        xpart_prior = F*xpart_prior + Gamma*gauss_rnd(0,sqrtQ^2,N);
        ypart = H * xpart_prior;
        for j = 1:N
            p_yx(1,j) = gauss_pdf(z_3(1, i),ypart(1,j),sqrtR^2).';
        end
        wgt = p_yx.*wgt;
        wgt = wgt/N;
        wgt = wgt/sum(wgt);

        xhat_SIS_3(:,i) = sum(xpart_prior.*repmat(wgt,2,1),2);
        P = (xpart_prior - repmat(xhat_SIS_3(:,i),1,N))*diag(wgt)*...
            (xpart_prior - repmat(xhat_SIS_3(:,i),1,N)).';
        sgmhat_SIS_3(1,i) = sqrt(P(1,1));
        sgmhat_SIS_3(2,i) = sqrt(P(2,2));
    end
    time_spent_SIS_3 = toc;
    
    % filtering (SIR)
    xhat_SIR_3 = zeros(2, simLen);
    sgmhat_SIR_3 = zeros(2, simLen);

    tic;
    wgt = 1/N*ones(1,N);
    xpart_prior = gauss_rnd(x0,sqrtP0^2,N);
    xhat_SIR_3(:,1) = sum(xpart_prior,2);
    sgmhat_SIR_3(:,1) = [sqrtP0(1,1);sqrtP0(2,2)];
    for i = 2:simLen
        resam_id = resampleSystematic(wgt);
        xpart_prior = xpart_prior(:,resam_id);
        xpart_prior = F*xpart_prior + Gamma*gauss_rnd(0,sqrtQ^2,N);

        ypart = H * xpart_prior;
        for j = 1:N
            p_yx(1,j) = gauss_pdf(z_3(1, i),ypart(1,j),sqrtR^2).';
        end
        wgt = p_yx.*wgt;
        wgt = wgt/N;
        wgt = wgt/sum(wgt);

        xhat_SIR_3(:,i) = sum(xpart_prior.*repmat(wgt,2,1),2);
        P = (xpart_prior - repmat(xhat_SIR_3(:,i),1,N))*diag(wgt)*...
            (xpart_prior - repmat(xhat_SIR_3(:,i),1,N)).';
        sgmhat_SIR_3(1,i) = sqrt(P(1,1));
        sgmhat_SIR_3(2,i) = sqrt(P(2,2));
    end
    time_spent_SIR_3 = toc;
    
% case 4 - sqrtQ = 0.01, sqrtR = 5
sqrtQ = 0.01;
sqrtR = 5;
    % simulation
    x_true_4 = zeros(2, simLen);
    z_true_4 = zeros(1, simLen);
    z_4 = zeros(1, simLen);
    x_true_4(:,1) = x0 + sqrtP0*randn(2,1);
    for i=2:simLen
        x_true_4(:,i) = F*x_true_4(:,i-1) + Gamma*sqrtQ*randn(1,1);
    end
    z_true_4(1,:) = x_true_4(1,:);
    z_4(1,:) = z_true_4 + sqrtR*randn(1,simLen);
    
    % filtering (SIS)
    xhat_SIS_4 = zeros(2, simLen);
    sgmhat_SIS_4 = zeros(2, simLen);

    tic;
    wgt = 1/N*ones(1,N);
    xpart_prior = gauss_rnd(x0,sqrtP0^2,N);
    xhat_SIS_4(:,1) = sum(xpart_prior,2);
    sgmhat_SIS_4(:,1) = [sqrtP0(1,1);sqrtP0(2,2)];
    for i = 2:simLen
        xpart_prior = F*xpart_prior + Gamma*gauss_rnd(0,sqrtQ^2,N);
        ypart = H * xpart_prior;
        for j = 1:N
            p_yx(1,j) = gauss_pdf(z_3(1, i),ypart(1,j),sqrtR^2).';
        end
        wgt = p_yx.*wgt;
        wgt = wgt/N;
        wgt = wgt/sum(wgt);

        xhat_SIS_4(:,i) = sum(xpart_prior.*repmat(wgt,2,1),2);
        P = (xpart_prior - repmat(xhat_SIS_4(:,i),1,N))*diag(wgt)*...
            (xpart_prior - repmat(xhat_SIS_4(:,i),1,N)).';
        sgmhat_SIS_4(1,i) = sqrt(P(1,1));
        sgmhat_SIS_4(2,i) = sqrt(P(2,2));
    end
    time_spent_SIS_4 = toc;
    
    % filtering (SIR)
    xhat_SIR_4 = zeros(2, simLen);
    sgmhat_SIR_4 = zeros(2, simLen);

    tic;
    wgt = 1/N*ones(1,N);
    xpart_prior = gauss_rnd(x0,sqrtP0^2,N);
    xhat_SIR_4(:,1) = sum(xpart_prior,2);
    sgmhat_SIR_4(:,1) = [sqrtP0(1,1);sqrtP0(2,2)];
    for i = 2:simLen
        resam_id = resampleSystematic(wgt);
        xpart_prior = xpart_prior(:,resam_id);
        xpart_prior = F*xpart_prior + Gamma*gauss_rnd(0,sqrtQ^2,N);

        ypart = H * xpart_prior;
        for j = 1:N
            p_yx(1,j) = gauss_pdf(z_4(1, i),ypart(1,j),sqrtR^2).';
        end
        wgt = p_yx.*wgt;
        wgt = wgt/N;
        wgt = wgt/sum(wgt);

        xhat_SIR_4(:,i) = sum(xpart_prior.*repmat(wgt,2,1),2);
        P = (xpart_prior - repmat(xhat_SIR_4(:,i),1,N))*diag(wgt)*...
            (xpart_prior - repmat(xhat_SIR_4(:,i),1,N)).';
        sgmhat_SIR_4(1,i) = sqrt(P(1,1));
        sgmhat_SIR_4(2,i) = sqrt(P(2,2));
    end
    time_spent_SIR_4 = toc;

%% results
% settings
default_font_size = 12.5; % default font size
default_line_width = 2; % default line width
     
% RMS error comparison
xerr_SIS_1 = x_true_1 - xhat_SIS_1;
xerr_SIS_2 = x_true_2 - xhat_SIS_2;
xerr_SIS_3 = x_true_3 - xhat_SIS_3;
xerr_SIS_4 = x_true_4 - xhat_SIS_4;
RMS_SIS_1 = [sqrt(sum(xerr_SIS_1(1,:).^2)/simLen);sqrt(sum(xerr_SIS_1(2,:).^2)/simLen)];
RMS_SIS_2 = [sqrt(sum(xerr_SIS_2(1,:).^2)/simLen);sqrt(sum(xerr_SIS_2(2,:).^2)/simLen)];
RMS_SIS_3 = [sqrt(sum(xerr_SIS_3(1,:).^2)/simLen);sqrt(sum(xerr_SIS_3(2,:).^2)/simLen)];
RMS_SIS_4 = [sqrt(sum(xerr_SIS_4(1,:).^2)/simLen);sqrt(sum(xerr_SIS_4(2,:).^2)/simLen)];
fprintf('RMSE (PF - SIS): \n')
fprintf('\t q = 10, r = 5: %f [m], %f [m/s]\n',RMS_SIS_1(1), RMS_SIS_1(2));
fprintf('\t q = 1, r = 5: %f [m], %f [m/s]\n',RMS_SIS_2(1), RMS_SIS_2(2));
fprintf('\t q = 0.1, r = 5: %f [m], %f [m/s]\n',RMS_SIS_3(1), RMS_SIS_3(2));
fprintf('\t q = 0.01, r = 5: %f [m], %f [m/s]\n',RMS_SIS_4(1), RMS_SIS_4(2));

xerr_SIR_1 = x_true_1 - xhat_SIR_1;
xerr_SIR_2 = x_true_2 - xhat_SIR_2;
xerr_SIR_3 = x_true_3 - xhat_SIR_3;
xerr_SIR_4 = x_true_4 - xhat_SIR_4;
RMS_SIR_1 = [sqrt(sum(xerr_SIR_1(1,:).^2)/simLen);sqrt(sum(xerr_SIR_1(2,:).^2)/simLen)];
RMS_SIR_2 = [sqrt(sum(xerr_SIR_2(1,:).^2)/simLen);sqrt(sum(xerr_SIR_2(2,:).^2)/simLen)];
RMS_SIR_3 = [sqrt(sum(xerr_SIR_3(1,:).^2)/simLen);sqrt(sum(xerr_SIR_3(2,:).^2)/simLen)];
RMS_SIR_4 = [sqrt(sum(xerr_SIR_4(1,:).^2)/simLen);sqrt(sum(xerr_SIR_4(2,:).^2)/simLen)];
fprintf('RMSE (PF - SIR): \n')
fprintf('\t q = 10, r = 5: %f [m], %f [m/s]\n',RMS_SIR_1(1), RMS_SIR_1(2));
fprintf('\t q = 1, r = 5: %f [m], %f [m/s]\n',RMS_SIR_2(1), RMS_SIR_2(2));
fprintf('\t q = 0.1, r = 5: %f [m], %f [m/s]\n',RMS_SIR_3(1), RMS_SIR_3(2));
fprintf('\t q = 0.01, r = 5: %f [m], %f [m/s]\n',RMS_SIR_4(1), RMS_SIR_4(2));
    
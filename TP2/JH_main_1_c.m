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
x_true = zeros(2, simLen);
z_true = zeros(1, simLen);
z = zeros(1, simLen);
x_true(:,1) = x0 + sqrtP0*randn(2,1);
for i=2:simLen
    x_true(:,i) = F*x_true(:,i-1) + Gamma*sqrtQ*randn(1,1);
end
z_true(1,:) = x_true(1,:);
z(1,:) = z_true + sqrtR*randn(1,simLen);

%% filtering
xtrue_array = reshape(repmat(reshape(x_true,2,[]),1,1),2,[]);

% PF - N = 100
xhat_SIS_100 = zeros(2, simLen);
sgmhat_SIS_100 = zeros(2, simLen);
N = 100; % particle numbers

tic;
wgt = 1/N*ones(1,N);
xpart_prior = gauss_rnd(x0,sqrtP0^2,N);
xhat_SIS_100(:,1) = sum(xpart_prior.*repmat(wgt,2,1),2);
sgmhat_SIS_100(:,1) = [sqrtP0(1,1);sqrtP0(2,2)];
for i = 2:simLen
    xpart_prior = F*xpart_prior + Gamma*gauss_rnd(0,sqrtQ^2,N);
    ypart = H * xpart_prior;
    for j = 1:N
        p_yx(1,j) = gauss_pdf(z(1, i),ypart(1,j),sqrtR^2).';
    end
    wgt = p_yx.*wgt;
    wgt = wgt/N;
    wgt = wgt/sum(wgt);

    xhat_SIS_100(:,i) = sum(xpart_prior.*repmat(wgt,2,1),2);
    P = (xpart_prior - repmat(xhat_SIS_100(:,i),1,N))*diag(wgt)*...
        (xpart_prior - repmat(xhat_SIS_100(:,i),1,N)).';
    sgmhat_SIS_100(1,i) = sqrt(P(1,1));
    sgmhat_SIS_100(2,i) = sqrt(P(2,2));
end
time_spent_SIS_100 = toc;
fprintf('100 Particles done ... \n');

% PF - N = 1000
xhat_SIS_1000 = zeros(2, simLen);
sgmhat_SIS_1000 = zeros(2, simLen);
N = 1000; % particle numbers

tic;
wgt = 1/N*ones(1,N);
xpart_prior = gauss_rnd(x0,sqrtP0^2,N);
xhat_SIS_1000(:,i) = sum(xpart_prior,2);
sgmhat_SIS_1000(:,i) = [sqrtP0(1,1);sqrtP0(2,2)];
for i = 2:simLen
    xpart_prior = F*xpart_prior + Gamma*gauss_rnd(0,sqrtQ^2,N);
    ypart = H * xpart_prior;
    for j = 1:N
        p_yx(1,j) = gauss_pdf(z(1, i),ypart(1,j),sqrtR^2).';
    end
    wgt = p_yx.*wgt;
    wgt = wgt/N;
    wgt = wgt/sum(wgt);

    xhat_SIS_1000(:,i) = sum(xpart_prior.*repmat(wgt,2,1),2);
    P = (xpart_prior - repmat(xhat_SIS_1000(:,i),1,N))*diag(wgt)*...
        (xpart_prior - repmat(xhat_SIS_1000(:,i),1,N)).';
    sgmhat_SIS_1000(1,i) = sqrt(P(1,1));
    sgmhat_SIS_1000(2,i) = sqrt(P(2,2));
end
time_spent_SIS_1000 = toc;
fprintf('1000 Particles done ... \n');

% PF - N = 10000
xhat_SIS_10000 = zeros(2, simLen);
sgmhat_SIS_10000 = zeros(2, simLen);
N = 10000; % particle numbers

tic;
wgt = 1/N*ones(1,N);
xpart_prior = gauss_rnd(x0,sqrtP0^2,N);
xhat_SIS_10000(:,i) = sum(xpart_prior,2);
sgmhat_SIS_10000(:,i) = [sqrtP0(1,1);sqrtP0(2,2)];
for i = 2:simLen
    xpart_prior = F*xpart_prior + Gamma*gauss_rnd(0,sqrtQ^2,N);
    ypart = H * xpart_prior;
    for j = 1:N
        p_yx(1,j) = gauss_pdf(z(1, i),ypart(1,j),sqrtR^2).';
    end
    wgt = p_yx.*wgt;
    wgt = wgt/N;
    wgt = wgt/sum(wgt);

    xhat_SIS_10000(:,i) = sum(xpart_prior.*repmat(wgt,2,1),2);
    P = (xpart_prior - repmat(xhat_SIS_10000(:,i),1,N))*diag(wgt)*...
        (xpart_prior - repmat(xhat_SIS_10000(:,i),1,N)).';
    sgmhat_SIS_10000(1,i) = sqrt(P(1,1));
    sgmhat_SIS_10000(2,i) = sqrt(P(2,2));
end
time_spent_SIS_10000 = toc;
fprintf('10000 Particles done ... \n\n');

%% results
% settings
default_font_size = 12.5; % default font size
default_line_width = 2; % default line width
     
% RMS error comparison
xerr_SIS_100 = xtrue_array - xhat_SIS_100;
xerr_SIS_1000 = xtrue_array - xhat_SIS_1000;
xerr_SIS_10000 = xtrue_array - xhat_SIS_10000;
RMS_SIS_100 = [sqrt(sum(xerr_SIS_100(1,:).^2)/simLen);sqrt(sum(xerr_SIS_100(2,:).^2)/simLen)];
RMS_SIS_1000 = [sqrt(sum(xerr_SIS_1000(1,:).^2)/simLen);sqrt(sum(xerr_SIS_1000(2,:).^2)/simLen)];
RMS_SIS_10000 = [sqrt(sum(xerr_SIS_10000(1,:).^2)/simLen);sqrt(sum(xerr_SIS_10000(2,:).^2)/simLen)];
fprintf('RMSE (PF - SIS): \n')
fprintf('\t N = 100: %f [m], %f [m/s]\n',RMS_SIS_100(1), RMS_SIS_100(2));
fprintf('\t N = 1000: %f [m], %f [m/s]\n',RMS_SIS_1000(1), RMS_SIS_1000(2));
fprintf('\t N = 10000: %f [m], %f [m/s]\n\n',RMS_SIS_10000(1), RMS_SIS_10000(2));
    
% computation time comparison
fprintf('Computation time (PF - SIS): \n')
fprintf('\t N = 100: %f [sec]\n',time_spent_SIS_100);
fprintf('\t N = 1000: %f [sec]\n',time_spent_SIS_1000);
fprintf('\t N = 10000: %f [sec]\n',time_spent_SIS_10000);

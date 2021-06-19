clear all; close all; clc

addpath('./resampling methods from Tiancheng Li');
addpath('./ekfukf');

%% Parameters
MC = 1;
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

% SIS 
xhat_SIS = zeros(MC*2, simLen);
sgmhat_SIS = zeros(MC*2, simLen);
time_spent_SIS = zeros(1, MC);
N = 10000; % particle numbers

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
    fprintf('%d-th MC done \n',mcid);
end
fprintf('PF - SIS done ... \n');

%% Data acquisition
save 1_c_10000 x_true xhat_SIS sgmhat_SIS time_spent_SIS
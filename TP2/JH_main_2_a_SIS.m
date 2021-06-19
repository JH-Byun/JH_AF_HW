%% Demonstration of univariate nonstationary growth model (UNGM) 
%
%  Description:
%    In this example various different non-linear filters and smoothers are
%    applied to the univariate nonstationary growth model (UNGM). The 
%    filters used in this demonstration are:
%      * Extended Kalman filter
%      * Unscented Kalman filter (augmented and non-augmented forms)
%      * Bootstrap filter
%      * Gauss-Hermite Kalman filter (degree 10)
%      * Cubature Kalman filter
%    Additionally, the corresponding Rauch-Tung-Striebel smoother results 
%    are also presented.
%
%  References:
%    Refer to the Toolbox documentation for details on the model.
%
%  See also:
%    ukf_predict1, ukf_update1, ukf_predict3, ukf_update3, urts_smooth1,
%    ekf_predict1, ekf_update1, erts_smooth1, ghkf_predict, ghkf_update, 
%    ghrts_smooth, ckf_predict, ckf_update, crts_smooth
%
%  Author:
%    Copyright (C) 2007 Jouni Hartikainen,
%    Updated by Arno Solin (2010).
%
%  Licence:
%    This software is distributed under the GNU General Public 
%    Licence (version 2 or later); please refer to the file 
%    Licence.txt, included with the software, for details.

%%
clear all; close all; clc

addpath('./resampling methods from Tiancheng Li');
addpath('./ekfukf');
addpath('./ekfukf/demos/ungm_demo')

%% Set up the model parameters

silent = 0;

clf; clc;
disp('This is a demonstration for EKF and UKF using the UNGM model.');
disp(' ');

% Handles to dynamic and measurement model functions,
% and their derivatives
f_func = @ungm_f;
h_func = @ungm_h;
df_func = @ungm_df_dx;
dh_func = @ungm_dh_dx;
d2f_func = @ungm_d2f_dx2;
d2h_func = @ungm_d2h_dx2;

% Number of samples
simLen = 500;

% number of Monte-Carlo samples
MC = 100;

% number of particle numbers
N = 500;

% Initial state and covariance
x_0 = .1;
P_0 = 2;

% Space for measurements
Y = zeros(1,simLen);

% Strengths of perturbations
u_n = 1;
v_n = 1;

fprintf('Generating real states and measurements...');
X = zeros(MC, simLen);
Y_r = zeros(MC, simLen);
for mcid = 1:MC
    X(mcid, 1) = ungm_f(x_0,1) + gauss_rnd(0,u_n,1);
    Y(mcid, 1) = feval(h_func,X(mcid, 1)) + gauss_rnd(0,v_n,1);
    for i = 2:simLen
        % Generate the true states with process noise
        X(mcid, i) = feval(f_func,X(mcid,i-1),i) + gauss_rnd(0,u_n,1);
        % Generate the observations with measurement noise
        Y(mcid, i) = feval(h_func,X(mcid, i)) + gauss_rnd(0,v_n,1);
    end
    Y_r(mcid,:) = feval(h_func,X(mcid,:));
end 
fprintf('Done!\n');

%% Run the various filters and smoothers
xhat = zeros(MC, simLen);
sgmhat = zeros(MC, simLen);
time_spent = zeros(1, MC);

for mcid = 1:MC
    tic;
    wgt = 1/N*ones(1,N);
    xpart_prior = gauss_rnd(x_0,P_0,N);
    xhat(mcid, 1) = sum(xpart_prior.*wgt,2);
    sgmhat(mcid,1) = sqrt((xpart_prior - repmat(xhat(mcid,1),1,N))*diag(wgt)*(xpart_prior - repmat(xhat(mcid,1),1,N)).');
    for i = 2:simLen
        xpart_prior = ungm_f(xpart_prior,i) + gauss_rnd(0,1,N);
        ypart = ungm_h(xpart_prior,i);
        for j = 1:N
            p_yx(1,j)=gauss_pdf(Y(mcid,i),ypart(1,j),1).'; % py_x
        end
        wgt = p_yx.*wgt;
        wgt = wgt/N;
        wgt = wgt/sum(wgt);
        
        xhat(mcid,i) = sum(xpart_prior.*wgt,2);
        sgmhat(mcid,i) = sqrt((xpart_prior - repmat(xhat(mcid,i),1,N))*diag(wgt)*(xpart_prior - repmat(xhat(mcid,i),1,N)).');
    end
    time_spent(1,mcid) = toc;
end

%% Data acquisition
save 2_a_SIS.mat X Y Y_r xhat sgmhat time_spent
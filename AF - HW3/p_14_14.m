clear all; close all; clc

%% simulation settings
t_f = 60;
T = 0.1;
tspan = 0:T:t_f;
F = [1 0 T 0;...
    0 1 0 T;...
    0 0 1 0;...
    0 0 0 1];
Q = diag([0 0 4 4]);
R = diag([1 1]);

% trajectory
x = zeros(4,length(tspan));
x(:,1) = [0 0 50 50].';
v(:,1) = mvnrnd(zeros(2,1),R);
y(:,1) = [sqrt((x(1,1)-20)^2+(x(2,1)-0)^2);...
        sqrt((x(1,1)-0)^2+(x(2,1)-20)^2)] + v(:,1);
for i = 1:length(tspan)-1
    w(:,i) = mvnrnd(zeros(4,1),Q);
    x(:,i+1) = F*x(:,i) + w(:,i);
    v(:,i+1) = mvnrnd(zeros(2,1),R);
    y(:,i+1) = [sqrt((x(1,i)-20)^2+(x(2,i)-0)^2);...
        sqrt((x(1,i)-0)^2+(x(2,i)-20)^2)] + v(:,i+1);
end

%% UKF
x_u = zeros(4,length(tspan));
x_u(:,1) = x(:,1);
P_u = zeros(4,4,length(tspan));
P_u(:,:,1) = zeros(4,4);

n = 4;
sgm_x = zeros(4,2*n);
sgm_y = zeros(2,2*n);
for i = 1:length(tspan) - 1
    % time update
    for j = 1:2*n
        temp = sqrtm(n*P_u(:,:,i));
        if j < n+1
            sgm_x(:,j) = x_u(:,i) + temp(j,:).';
        else
            sgm_x(:,j) = x_u(:,i) - temp(j-n,:).';
        end
        sgm_x(:,j) = F*sgm_x(:,j);
    end
    
    x_priori = zeros(4,1);
    for j = 1:2*n
        x_priori = x_priori + (1/(2*n))*sgm_x(:,j);
    end
    
    P_priori = zeros(4,4);
    for j = 1:2*n
        P_priori = P_priori + ...
            (1/(2*n))*((sgm_x(:,j)-x_priori)*(sgm_x(:,j)-x_priori).'+Q);
    end
    
    % measurement update
    for j = 1:2*n
        temp = sqrtm(n*P_priori);
        if j < n+1
            sgm_x(:,j) = x_priori + temp(j,:).';
        else
            sgm_x(:,j) = x_priori - temp(j-n,:).';
        end
        sgm_y(:,j) = [sqrt((sgm_x(1,j)-20)^2+(sgm_x(2,j)-0)^2);...
            sqrt((sgm_x(1,j)-0)^2+(sgm_x(2,j)-20)^2)];
    end
    
    y_u = zeros(2,1);
    for j = 1:2*n
        y_u = y_u + (1/(2*n))*sgm_y(:,j);
    end
    
    P_y = zeros(2,2);
    P_xy = zeros(4,2);
    for j = 1:2*n
        P_y = P_y + (1/(2*n))*((sgm_y(:,j)-y_u)*(sgm_y(:,j)-y_u).'+R);
        P_xy = P_xy + (1/(2*n))*(sgm_x(:,j)-x_priori)*(sgm_y(:,j)-y_u).';
    end
    K = P_xy*inv(P_y);
    x_u(:,i+1) = x_priori + K*(y(:,i+1)-y_u);
    P_u(:,:,i+1) = P_priori - K*P_y*K.';
end

% plot (UKF)
linewidth = 2;
figure;
subplot(2,2,1)
plot(tspan,x_u(1,:),'LineWidth',linewidth)
hold on
plot(tspan,x(1,:),'LineWidth',linewidth)
grid on
xlabel('time [sec]');
ylabel('n [m]');
subplot(2,2,2)
plot(tspan,x_u(2,:),'LineWidth',linewidth)
hold on
plot(tspan,x(2,:),'LineWidth',linewidth)
grid on
xlabel('time [sec]');
ylabel('e [m]');
subplot(2,2,3)
plot(tspan,x_u(3,:),'LineWidth',linewidth)
hold on
plot(tspan,x(3,:),'LineWidth',linewidth)
grid on
xlabel('time [sec]');
ylabel('$\dot{n}$ [m/s]','interpreter','latex');
subplot(2,2,4)
plot(tspan,x_u(4,:),'LineWidth',linewidth)
hold on
plot(tspan,x(4,:),'LineWidth',linewidth)
grid on
xlabel('time [sec]');
ylabel('$\dot{e}$ [m/s]','interpreter','latex');
legend('estimated (ukf)','true');

figure;
subplot(2,2,1)
plot(tspan,x_u(1,:)-x(1,:),'LineWidth',linewidth)
grid on
xlabel('time [sec]');
ylabel('$\hat{n}-n$ [m]','interpreter','latex');
subplot(2,2,2)
plot(tspan,x_u(2,:)-x(2,:),'LineWidth',linewidth)
grid on
xlabel('time [sec]');
ylabel('$\hat{e}-e$ [m]','interpreter','latex');
subplot(2,2,3)
plot(tspan,x_u(3,:)-x(3,:),'LineWidth',linewidth)
grid on
xlabel('time [sec]');
ylabel('$\dot{\hat{n}}-\dot{n}$ [m/s]','interpreter','latex');
subplot(2,2,4)
plot(tspan,x_u(4,:)-x(4,:),'LineWidth',linewidth)
grid on
xlabel('time [sec]');
ylabel('$\dot{\hat{e}}-\dot{e}$ [m/s]','interpreter','latex');

% RMSE (ekf)
MSE = zeros(4,1);
for i = 1:length(tspan)
    MSE = MSE + (x_u(:,i)-x(:,i)).^2;
end
RMSE_u(1) = sqrt(MSE(1)/length(tspan));
RMSE_u(2) = sqrt(MSE(2)/length(tspan));
RMSE_u(3) = sqrt(MSE(3)/length(tspan));
RMSE_u(4) = sqrt(MSE(4)/length(tspan));

%% EKF
x_ekf = zeros(4,length(tspan));
x_ekf(:,1) = x(:,1);
P_ekf = zeros(4,4,length(tspan));
P_ekf(:,:,1) = zeros(4,4);

for i = 1:length(tspan)-1
    % time update
    P_priori = F*P_ekf(:,:,i)*F.' + Q;
    x_priori = F*x_ekf(:,i);
    % measurement update
    n = x_priori(1);
    e = x_priori(2);
    H = zeros(2,4);
    H(1,1) = (n-20)/sqrt((n-20)^2+e^2);
    H(1,2) = e/sqrt((n-20)^2+e^2);
    H(2,1) = n/sqrt(n^2+(e-20)^2);
    H(2,2) = (e-20)/sqrt(n^2+(e-20)^2);
    K = P_priori*H.'*inv(H*P_priori*H.'+R);
    x_ekf(:,i+1) = x_priori + K*(y(:,i+1)-[sqrt((n-20)^2+e^2);sqrt(n^2+(e-20)^2)]);
    P_ekf(:,:,i+1) = (eye(4)-K*H)*P_priori;
end

% plot (EKF)
linewidth = 2;
figure;
subplot(2,2,1)
plot(tspan,x_ekf(1,:),'LineWidth',linewidth)
hold on
plot(tspan,x(1,:),'LineWidth',linewidth)
grid on
xlabel('time [sec]');
ylabel('n [m]');
subplot(2,2,2)
plot(tspan,x_ekf(2,:),'LineWidth',linewidth)
hold on
plot(tspan,x(2,:),'LineWidth',linewidth)
grid on
xlabel('time [sec]');
ylabel('e [m]');
subplot(2,2,3)
plot(tspan,x_ekf(3,:),'LineWidth',linewidth)
hold on
plot(tspan,x(3,:),'LineWidth',linewidth)
grid on
xlabel('time [sec]');
ylabel('$\dot{n}$ [m/s]','interpreter','latex');
subplot(2,2,4)
plot(tspan,x_ekf(4,:),'LineWidth',linewidth)
hold on
plot(tspan,x(4,:),'LineWidth',linewidth)
grid on
xlabel('time [sec]');
ylabel('$\dot{e}$ [m/s]','interpreter','latex');
legend('estimated (ekf)','true');

figure;
subplot(2,2,1)
plot(tspan,x_ekf(1,:)-x(1,:),'LineWidth',linewidth)
grid on
xlabel('time [sec]');
ylabel('$\hat{n}-n$ [m]','interpreter','latex');
subplot(2,2,2)
plot(tspan,x_ekf(2,:)-x(2,:),'LineWidth',linewidth)
grid on
xlabel('time [sec]');
ylabel('$\hat{e}-e$ [m]','interpreter','latex');
subplot(2,2,3)
plot(tspan,x_ekf(3,:)-x(3,:),'LineWidth',linewidth)
grid on
xlabel('time [sec]');
ylabel('$\dot{\hat{n}}-\dot{n}$ [m/s]','interpreter','latex');
subplot(2,2,4)
plot(tspan,x_ekf(4,:)-x(4,:),'LineWidth',linewidth)
grid on
xlabel('time [sec]');
ylabel('$\dot{\hat{e}}-\dot{e}$ [m/s]','interpreter','latex');

% RMSE (ekf)
MSE = zeros(4,1);
for i = 1:length(tspan)
    MSE = MSE + (x_ekf(:,i)-x(:,i)).^2;
end
RMSE_ekf(1) = sqrt(MSE(1)/length(tspan));
RMSE_ekf(2) = sqrt(MSE(2)/length(tspan));
RMSE_ekf(3) = sqrt(MSE(3)/length(tspan));
RMSE_ekf(4) = sqrt(MSE(4)/length(tspan));
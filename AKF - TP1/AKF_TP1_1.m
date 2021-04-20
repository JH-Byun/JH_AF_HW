clear all; close all; clc

%% simulation settings
N = 200;
dt = 0.1;
x = zeros(6, N+1);
x(1,:) = normnd([0 0 0 -1 0 0].',diag([0.1^2 0.1^2 0.1^2 0.1^2 0.5^2 0.5^2]));

F_1 = [0 0 1 0;...
    0 0 0 1;...
    0 0 0 0;...
    0 0 0 0];
L_1 = [zeros(2,2);eye(2)];
Phi_1 = eye(4) + F_1*dt;

F_2 = [zeros(4,2) eye(4);zeros(2,6)];
L_2 = [zeros(2,2);eye(2)];
Phi_2 = eye(6) + F_2*dt + F_2^2*(dt^2/2);

%%
A = eye(4) + F*dt;
Q = (1/2)*L*Q_*L.';


for i = 1:N
    
end



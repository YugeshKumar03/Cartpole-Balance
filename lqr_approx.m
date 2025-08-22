clc, clearvars, close all

% Defining constants
m = 0.1;
M = 1;
g = 9.81;
l = 1;

A_m = zeros(4,4);
A_m(1,2) = 1;
A_m(2,3) = m*g/M;
A_m(3,4) = 1;
A_m(4,3) = g*(1 + m/M)/l;

B_m = zeros(4, 1);
B_m(2,1) = 1/M;
B_m(4,1) = 1/(M*l);

if rank(ctrb(A_m, B_m)) == 4
    fprintf("This System is controllable!\n");
else
    fprintf("This System is not controllable..n");
end

% Q_m = eye(4);
Q_m = [1 0 0 0;
       0 1 0 0;
       0 0 1 0;
       0 0 0 1];
R = 0.01;

n = 10;
t_sim = linspace(0, n, 10*n+1)'; % Simulation Time

x_i = [0; 0; 20*pi/180; 0]; % Initial conditions on state vector
x_ref = [0; 0; 0; 0]; % Equilibrium position

% Finding the matrix K_m which minimizes the defined weight function
K_m = lqr(A_m, B_m, Q_m, R);

% Equivalent system of form X^. = A1_m*X
A1_m = A_m-B_m*K_m;
[T_m, D_m] = eig(A1_m);

X_f = zeros(length(t_sim), 4);
% X_f(1, :) = x_i;

for i = 1:length(t_sim)
    %X_f(i, :) = (T_m*exp(D_m*t_sim(i))*(T_m\(x_i-x_ref)) + x_ref)';
    X_f(i, :) = (x_ref + T_m*expm(t_sim(i)*D_m)*(T_m\(x_i-x_ref)))';
end

% X_f = T_m*exp(D_m*t_sim)*(T_m\x_i);
% plot(t_sim, X_f);

X_f = real(X_f);

% figure(2)
% plot(t_sim, X_f(:,3));
% plot(t_sim, -K_m*X_f');

CtP_render(t_sim, X_f, m, M, l);


%% ASEN 5014 - Linear Control Systems
% Final Project
% Galen Savidge, Aniket Goel, Andrew Palski

clear; close all; format shortG; clc;

% Question 2:
% Linear system
[A, B, C, D, G, B_tot, D_tot] = sys_setup();

% Initial condition
%x0 = [0; 10; 0; 0; 0; .001]; 
%x0 = [-.75; 5; 0; -.001; 0; .001];
x0 = [0; 5; 0; -.001; 0; .001];

r = [0; 0.5; 0]; % Reference input [km]
d = -1e-9; % Disturbance of 1 um/s^2 [km/s^2]
umax = 10 / 1300 * 1e-3; % Maximum acceleration per thruster [km/s^2]

% Sim setup
ts = 0:1:18000;
rs = repmat(r, 1, length(ts)); % Reference input history
us = repmat([r', d],length(ts),1); % Reference input augmented with disturbance
[K, F] = fsf_gains(A, B, C);

% Question 5
% Luenberger observer matrix
des_L_poles = [-3, -3.1, -3.2, -3.3, -3.4, -3.5]*0.01;
L_transpose = place(A',C',des_L_poles);
L = L_transpose';
% Closed loop state and error dynamics with full state feedback control
A_CL_aug = [(A-B*K), B*K;
        zeros(size(A-B*K)), (A-L*C)];
B_CL_aug = [B*F; zeros(size(B*F))];
C_CL_aug = [C, zeros(size(C))];
D_CL_aug = D_tot;
B_CL_aug_tot = [B_CL_aug, [G;G]];

% Simulating with zero initial error
e0 = zeros(6,1);
x_cl_aug0 = [x0;e0];

sys_CL_aug = ss(A_CL_aug,B_CL_aug_tot,C_CL_aug,D_CL_aug);
[ys,~,xs] = lsim(sys_CL_aug,us,ts,x_cl_aug0);

plot_state(ts, xs, 'Simulated States (With Observer, e(0) = 0)')
plot_state(ts, xs(:,7:12), 'Observer Error (e(0) = 0)')

% Simulating with nonzero initial error
e0 = [1 1 1 .1 .1 .1]';
x_cl_aug0 = [x0;e0];

sys_CL_aug = ss(A_CL_aug,B_CL_aug_tot,C_CL_aug,D_CL_aug);
[ys,~,xs] = lsim(sys_CL_aug,us,ts,x_cl_aug0);

plot_state(ts, xs, 'Simulated States (With Observer, e(0) \neq 0)')
plot_state(ts, xs(:,7:12), 'Observer Error (e(0) \neq 0)')

% Question 6
% Set up cost function using Bryson's rules
xmax = ones(1, 6); % TBD

% Tuning parameters
alpha = ones(1, 6); % State error weights
alpha = alpha./sum(alpha);
beta = ones(1, 3); % Input weights
beta = beta./sum(beta);
rho = 2e7;

Q = diag(alpha./(xmax.^2))
R = rho*diag(beta/umax)

% Use LQR to generate CL gain from cost function
[K_LQR,W,evals_LQR] = lqr(A,B,Q,R);

% Feed-forward gain for zero steady-state offset
F_LQR = inv(C/(-A+B*K_LQR)*B);

% Assemble closed-loop system (with disturbance and observer)
Acl_LQR = [A-B*K_LQR, B*K_LQR;
           zeros(size(A-B*K_LQR)), A-L*C];
Bcl_LQR = [B*F_LQR, G;
           zeros(size(B*F_LQR)), G];
Ccl_LQR = [C, zeros(size(C))];
Dcl_LQR = D_tot;
syscl_LQR = ss(Acl_LQR, Bcl_LQR, Ccl_LQR, Dcl_LQR);

% Simulate system
[ys,~,xs] = lsim(syscl_LQR,us,ts,[x0; zeros(6,1)]);

plot_state(ts, xs, 'Simualted States (LQR)')
plot_actuator_responses(ts, rs, xs, F_LQR, K_LQR, umax, 'Actuator Responses (LQR)')

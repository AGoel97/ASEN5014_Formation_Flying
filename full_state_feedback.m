%% ASEN 5014 - Linear Control Systems
% Final Project Part B Question 4: Full State Feedback
% Galen Savidge, Aniket Goel, Andrew Palski

clear; close all; format shortG; clc;

% Linear system
[A, B, C, D, G, B_tot, D_tot] = sys_setup();

% Initial condition
x0 = [0; 5; 0; -.001; 0; .001];

r = [0; 0.5; 0]; % Reference input [km]
d = -1e-9; % Disturbance of 1 um/s^2 [km/s^2]
umax = 10 / 1300 * 1e-3; % Maximum acceleration per thruster [km/s^2]

% Sim setup
ts = 0:1:18000;
rs = repmat(r, 1, length(ts)); % Reference input history

% Try without integral control at first
[K, F] = fsf_gains(A, B, C);

Acl_FSF = A - B*K;
Bcl_FSF = [B*F G];
Ccl_FSF = C - D*K;
Dcl_FSF = [D*F zeros(3,1)];

sys_FSF = ss(Acl_FSF, Bcl_FSF, Ccl_FSF, Dcl_FSF);

% Open-loop augmented system for integral control
Aaug = [A zeros(6,3);-C zeros(3)];
Baug = [B; zeros(3)];
Caug = [C zeros(3)];
Daug = zeros(3);
Faug = [zeros(6,3);eye(3)];
Gaug = [G; zeros(3,1)];

r_aug_FSF = repmat([r' d],length(ts),1);%reference input augmented with disturbance

[y_FSF,~,x_FSF] = lsim(sys_FSF, r_aug_FSF, ts, x0);

plot_state(ts, x_FSF, 'Simulated States (Full State Feedback)')
plot_actuator_responses(ts, rs, x_FSF, F, K, umax, 'Actuator Responses (Full State Feedback)')

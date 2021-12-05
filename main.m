%% ASEN 5014 - Linear Control Systems
% Final Project
% Galen Savidge, Aniket Goel, Andrew Palski

clear; close all; format shortG; clc;
%% PART A:Linear Systems Analysis
% Question 2:
% Linear system
n = sqrt(398600/6778^3);
A = [0 0 0 1 0 0;
     0 0 0 0 1 0;
     0 0 0 0 0 1;
     3*n^2 0 0 0 2*n 0;
     0 0 0 -2*n 0 0;
     0 0 -n^2 0 0 0];
B = [0 0 0;
     0 0 0;
     0 0 0;
     1 0 0;
     0 1 0;
     0 0 1];
C = [1 0 0 0 0 0;
     0 1 0 0 0 0;
     0 0 1 0 0 0];
D = [0 0 0;
     0 0 0;
     0 0 0];
G = [0;0;0;0;1;0]; %disturbance matrix, disturbance only affects in track dynamics
B_tot = [B,G];
D_tot = [D,[0;0;0]];
x0 = [0; 10; 0; 0; 0; .001]; % initial condition
r = [0; 5; 0; 0; 0; 0]; % reference 
d = -1e-9; % disturbance %CHANGE: set to 1 um/s vs 1 m/s

 % Eigenvalues
[E, L] = eig(A);
evals = diag(L);

% Simulated plant dynamics with no control and including disturbances
u = [0;0;0];
u_tot = [u;d]';
ts = 0:1:18000;
us = u_tot + zeros(length(ts),4);
sys = ss(A,B_tot,C,D_tot);
[ys,~,xs] = lsim(sys,us,ts,x0);
figure();
ax = subplot(3,1,1);
plot(ax,ts,xs(:,1),'LineWidth',2,'Color','r')
ylabel('x - radial position (km)')
title('Simulated States (Positions)')
grid on
ax = subplot(3,1,2);
plot(ax,ts,xs(:,2),'LineWidth',2,'Color','k')
ylabel('y - in-track position (km)')
grid on
ax = subplot(3,1,3);
plot(ax,ts,xs(:,3),'LineWidth',2,'Color','b')
ylabel('z - cross-track position (km)')
xlabel('Time (sec)')
grid on

figure()
ax = subplot(3,1,1);
plot(ax,ts,xs(:,4),'LineWidth',2,'Color','r')
ylabel('xdot - radial velocity (km/s)')
title('Simulated States (Velocities)')
grid on
ax = subplot(3,1,2);
plot(ax,ts,xs(:,5),'LineWidth',2,'Color','k')
ylabel('ydot - in-track velocity (km/s)')
grid on
ax = subplot(3,1,3);
plot(ax,ts,xs(:,6),'LineWidth',2,'Color','b')
ylabel('zdot - cross-track velocity (km/s)')
xlabel('Time (sec)')
grid on

% Question 3
% Reachability/Controllability
P = [B A*B A^2*B A^3*B A^4*B A^5*B]; 
rank(P) % = 6 = n, so reachable

% Observability
O = [C; C*A; C*A^2; C*A^3; C*A^4; C*A^5];
rank(O) % = 6 = n, so completely observable


% Question 5
% Luenberger observer matrix
des_L_poles = [-2, -2.1, -2.2, -2.3, -2.4, -2.5];
L_transpose = place(A',C',des_L_poles);
% Closed loop state and error dynamics with full state feedback control
A_CL_aug = [(A-B*K), B*K;
        0, (A-L*C)];
B_CL_aug = [B*F; zeros(size(B*F))];
C_CL_aug = [C, zeros(size(C))];
D_CL_aug = D_tot;
B_CL_aug_tot = [B_CL_aug, [G;G]];

% Simulating with zero initial error
e0 = zeros(6,1);
x_cl_aug0 = [x0;e0];

u_tot = [r;d];
ts = 0:1:18000;
us = u_tot + zeros(length(ts),4);
sys = ss(A_CL_aug,B_CL_aug_tot,C_CL_aug,D_CL_aug);
[ys,~,xs] = lsim(sys,us,ts,x0);
figure();
ax = subplot(3,1,1);
plot(ax,ts,xs(:,1),'LineWidth',2,'Color','r')
ylabel('x - radial position (km)')
title('Simulated States (Positions)')
grid on
ax = subplot(3,1,2);
plot(ax,ts,xs(:,2),'LineWidth',2,'Color','k')
ylabel('y - in-track position (km)')
grid on
ax = subplot(3,1,3);
plot(ax,ts,xs(:,3),'LineWidth',2,'Color','b')
ylabel('z - cross-track position (km)')
xlabel('Time (sec)')
grid on

figure()
ax = subplot(3,1,1);
plot(ax,ts,xs(:,4),'LineWidth',2,'Color','r')
ylabel('xdot - radial velocity (km/s)')
title('Simulated States (Velocities)')
grid on
ax = subplot(3,1,2);
plot(ax,ts,xs(:,5),'LineWidth',2,'Color','k')
ylabel('ydot - in-track velocity (km/s)')
grid on
ax = subplot(3,1,3);
plot(ax,ts,xs(:,6),'LineWidth',2,'Color','b')
ylabel('zdot - cross-track velocity (km/s)')
xlabel('Time (sec)')
grid on

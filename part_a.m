%% ASEN 5014 - Linear Control Systems
% Final Project Part A: Linear Systems Analysis
% Galen Savidge, Aniket Goel, Andrew Palski

clear; close all; format shortG; clc;

% Question 2:
% Linear system
[A, B, C, D, G, B_tot, D_tot] = sys_setup();

% Initial condition
%x0 = [0; 10; 0; 0; 0; .001]; 
%x0 = [-.75; 5; 0; -.001; 0; .001];
x0 = [0; 5; 0; -.001; 0; .001];

d = -1e-9; % Disturbance of 1 um/s^2 [km/s^2]

% Eigenvalues
[E, L] = eig(A);
evals = diag(L);

% Simulated plant dynamics with no control and including disturbances
u = [0;0;0];
u_tot = [u;d]';
ts = 0:1:18000;
us = u_tot + zeros(length(ts),4);
sys_OL = ss(A,B_tot,C,D_tot);
[ys,~,xs] = lsim(sys_OL,us,ts,x0);

% Plot open-loop system response
figure()
fig = gcf;
fig.Position = [0 50 1000 650];
ax = subplot(3,2,1);
plot(ax,ts,xs(:,1),'LineWidth',2,'Color','r')
ylabel('x - radial position (km)')
grid on
ax = subplot(3,2,3);
plot(ax,ts,xs(:,2),'LineWidth',2,'Color','k')
ylabel('y - in-track position (km)')
grid on
ax = subplot(3,2,5);
plot(ax,ts,xs(:,3),'LineWidth',2,'Color','b')
ylabel('z - cross-track position (km)')
xlabel('Time (sec)')
grid on
ax = subplot(3,2,2);
plot(ax,ts,xs(:,4),'LineWidth',2,'Color','r')
ylabel('xdot - radial velocity (km/s)')
grid on
ax = subplot(3,2,4);
plot(ax,ts,xs(:,5),'LineWidth',2,'Color','k')
ylabel('ydot - in-track velocity (km/s)')
grid on
ax = subplot(3,2,6);
plot(ax,ts,xs(:,6),'LineWidth',2,'Color','b')
ylabel('zdot - cross-track velocity (km/s)')
xlabel('Time (sec)')
grid on
sgtitle('Open-loop System Response')

% Question 3
% Reachability/Controllability
P = [B A*B A^2*B A^3*B A^4*B A^5*B]; 
rank(P) % = 6 = n, so reachable

% Observability
O = [C; C*A; C*A^2; C*A^3; C*A^4; C*A^5];
rank(O) % = 6 = n, so completely observable
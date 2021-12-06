%% ASEN 5014 - Linear Control Systems
% Final Project
% Galen Savidge, Aniket Goel, Andrew Palski

clear; close all; format shortG; clc;
%% PART A: Linear Systems Analysis
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
sys_OL = ss(A,B_tot,C,D_tot);
[ys,~,xs] = lsim(sys_OL,us,ts,x0);
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

%% PART B
% Question 4: Full State Feedback

% Try without integral control at first

K = place(A,B,[-0.0013149 -0.00327 -0.00156 -.00301 -.00264 -.000828]);
%K = place(A,B,-.005*rand(1,6));
F = inv(C/(-A+B*K)*B);

Acl_FSF = A - B*K;
Bcl_FSF = [B*F G];
Ccl_FSF = C - D*K;
Dcl_FSF = [D*F zeros(3,1)];

sys_FSF = ss(Acl_FSF, Bcl_FSF, Ccl_FSF, Dcl_FSF);


r_aug_FSF = repmat([0 .5 0 d],length(ts),1);%reference input augmented with disturbance

[y_FSF,~,x_FSF] = lsim(sys_FSF, r_aug_FSF, ts, x0);

%Check actuators
u = zeros(3,length(ts));
for ii = 1:length(ts)
    u(:,ii) = F*[0;.5;0] - K * x_FSF(ii,:)';
end

u = 1000 * u; %convert to m/s

figure();
ax = subplot(3,1,1);
plot(ax,ts,x_FSF(:,1),'LineWidth',2,'Color','r')
ylabel('x - radial position (km)')
title('Simulated States (Positions)')
grid on
ax = subplot(3,1,2);
plot(ax,ts,x_FSF(:,2),'LineWidth',2,'Color','k')
ylabel('y - in-track position (km)')
grid on
ax = subplot(3,1,3);
plot(ax,ts,x_FSF(:,3),'LineWidth',2,'Color','b')
ylabel('z - cross-track position (km)')
xlabel('Time (sec)')
grid on

figure()
ax = subplot(3,1,1);
plot(ax,ts,x_FSF(:,4),'LineWidth',2,'Color','r')
ylabel('xdot - radial velocity (km/s)')
title('Simulated States (Velocities)')
grid on
ax = subplot(3,1,2);
plot(ax,ts,x_FSF(:,5),'LineWidth',2,'Color','k')
ylabel('ydot - in-track velocity (km/s)')
grid on
ax = subplot(3,1,3);
plot(ax,ts,x_FSF(:,6),'LineWidth',2,'Color','b')
ylabel('zdot - cross-track velocity (km/s)')
xlabel('Time (sec)')
grid on

figure()
subplot(3,1,1);
plot(ts,u(1,:),'LineWidth',2,'Color','r');
hold on;
plot([0 max(ts)],[0.000769 0.000769],'k:');
plot([0 max(ts)],[-0.000769 -0.000769],'k:');
ylabel('Radial Actuator Response (m/s^2)');

subplot(3,1,2);
plot(ts,u(2,:),'LineWidth',2,'Color','k');
hold on;
plot([0 max(ts)],[0.000769 0.000769],'k:');
plot([0 max(ts)],[-0.000769 -0.000769],'k:');
ylabel('In-Track Actuator Response (m/s^2)');

subplot(3,1,3);
plot(ts,u(3,:),'LineWidth',2,'Color','b');
hold on;
plot([0 max(ts)],[0.000769 0.000769],'k:');
plot([0 max(ts)],[-0.000769 -0.000769],'k:');
ylabel('Cross-Track Actuator Response (m/s^2)');
xlabel('Time (sec)');
sgtitle('Actuator Response vs Time');

% Question 5
% Luenberger observer matrix
des_L_poles = [-3, -3.1, -3.2, -3.3, -3.4, -3.5];
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

us = r_aug_FSF;
sys_CL_aug = ss(A_CL_aug,B_CL_aug_tot,C_CL_aug,D_CL_aug);
[ys,~,xs] = lsim(sys_CL_aug,us,ts,x_cl_aug0);
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

figure();
ax = subplot(3,1,1);
plot(ax,ts,xs(:,7),'LineWidth',2,'Color','r')
ylabel('x - radial position (km)')
title('Error (Positions)')
grid on
ax = subplot(3,1,2);
plot(ax,ts,xs(:,8),'LineWidth',2,'Color','k')
ylabel('y - in-track position (km)')
grid on
ax = subplot(3,1,3);
plot(ax,ts,xs(:,9),'LineWidth',2,'Color','b')
ylabel('z - cross-track position (km)')
xlabel('Time (sec)')
grid on

figure();
ax = subplot(3,1,1);
plot(ax,ts,xs(:,10),'LineWidth',2,'Color','r')
ylabel('x - radial position (km)')
title('Error (Velocities)')
grid on
ax = subplot(3,1,2);
plot(ax,ts,xs(:,11),'LineWidth',2,'Color','k')
ylabel('y - in-track position (km)')
grid on
ax = subplot(3,1,3);
plot(ax,ts,xs(:,12),'LineWidth',2,'Color','b')
ylabel('z - cross-track position (km)')
xlabel('Time (sec)')
grid on

% Simulating with nonzero initial error
e0 = [1 1 1 .1 .1 .1]';
x_cl_aug0 = [x0;e0];

us = r_aug_FSF;
sys_CL_aug = ss(A_CL_aug,B_CL_aug_tot,C_CL_aug,D_CL_aug);
[ys,~,xs] = lsim(sys_CL_aug,us,ts,x_cl_aug0);
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

figure();
ax = subplot(3,1,1);
plot(ax,ts,xs(:,7),'LineWidth',2,'Color','r')
ylabel('x - radial position (km)')
title('Error (Positions)')
grid on
ax = subplot(3,1,2);
plot(ax,ts,xs(:,8),'LineWidth',2,'Color','k')
ylabel('y - in-track position (km)')
grid on
ax = subplot(3,1,3);
plot(ax,ts,xs(:,9),'LineWidth',2,'Color','b')
ylabel('z - cross-track position (km)')
xlabel('Time (sec)')
grid on

figure();
ax = subplot(3,1,1);
plot(ax,ts,xs(:,10),'LineWidth',2,'Color','r')
ylabel('x - radial position (km)')
title('Error (Velocities)')
grid on
ax = subplot(3,1,2);
plot(ax,ts,xs(:,11),'LineWidth',2,'Color','k')
ylabel('y - in-track position (km)')
grid on
ax = subplot(3,1,3);
plot(ax,ts,xs(:,12),'LineWidth',2,'Color','b')
ylabel('z - cross-track position (km)')
xlabel('Time (sec)')

% Question 6
% Set up cost function using Bryson's rules
xmax = ones(1, 6); % TBD
umax = 1 / 1300 * 1e-3 % Maximum acceleration per thruster [km/s^2]

% Tuning parameters
alpha = ones(1, 6); % State error weights
alpha = alpha./sum(alpha);
beta = ones(1, 3); % Input weights
beta = beta./sum(beta);
rho = 2e6;

Q = diag(alpha./(xmax.^2))
R = rho*diag(beta/umax)

% Use LQR to generate CL gain from cost function
[K_LQR,W,evals_LQR] = lqr(A,B,Q,R);

% Feed-forward gain for zero steady-state offset
F_LQR = inv(C/(-A+B*K_LQR)*B);

% Assemble closed-loop system (with disturbance)
Acl_LQR = A - B*K_LQR;
Bcl_LQR = [B*F_LQR, G];
Ccl_LQR = C;
Dcl_LQR = [D, zeros(3,1)];
syscl_LQR = ss(Acl_LQR, Bcl_LQR, Ccl_LQR, Dcl_LQR);

% Simulate system
[ys,~,xs] = lsim(syscl_LQR,r_aug_FSF,ts,x0);

% Reconstruct inputs
rs = repmat([0 0.5 0], length(ts), 1);
u = F_LQR*rs' - K_LQR*xs';
u = u * 1e3;

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

figure()
subplot(3,1,1);
plot(ts,u(1,:),'LineWidth',2,'Color','r');
hold on;
plot([0 max(ts)],[0.000769 0.000769],'k:');
plot([0 max(ts)],[-0.000769 -0.000769],'k:');
ylabel('Radial Actuator Response (m/s^2)');

subplot(3,1,2);
plot(ts,u(2,:),'LineWidth',2,'Color','k');
hold on;
plot([0 max(ts)],[0.000769 0.000769],'k:');
plot([0 max(ts)],[-0.000769 -0.000769],'k:');
ylabel('In-Track Actuator Response (m/s^2)');

subplot(3,1,3);
plot(ts,u(3,:),'LineWidth',2,'Color','b');
hold on;
plot([0 max(ts)],[0.000769 0.000769],'k:');
plot([0 max(ts)],[-0.000769 -0.000769],'k:');
ylabel('Cross-Track Actuator Response (m/s^2)');
xlabel('Time (sec)');
sgtitle('Actuator Response vs Time');
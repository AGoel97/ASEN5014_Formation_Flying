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

%% Part B
% Question 4: Full State Feedback

% Try without integral control at first

%K = place(A,B,...
    %[-0.0010652 -0.0016794 -0.0013118 -.0018838 -.00088127 -.00090389]);
K = place(A,B,-.002*rand(1,6));
F = inv(C/(-A+B*K)*B);

Acl_FSF = A - B*K;
Bcl_FSF = [B*F G];
Ccl_FSF = C - D*K;
Dcl_FSF = [D*F zeros(3,1)];

sys_FSF = ss(Acl_FSF, Bcl_FSF, Ccl_FSF, Dcl_FSF);

% Having issues getting integral control to not overshoot
%With integral control, first set up augmented matrices

% Aaug_FSF = [A zeros(6,3);-C zeros(3)];
% Baug_FSF = [B; zeros(3)];
% Caug_FSF = [C zeros(3)];
% Daug_FSF = zeros(3);
% Faug_FSF = [zeros(6,3);eye(3)];
% Gaug_FSF = [G; zeros(3,1)];
% 
% Kaug_FSF = place(Aaug_FSF,Baug_FSF,...
%     -.003*rand(1,9));
%     %[-0.0013149 -0.00327 -0.00156 -.00301 -.00264 -.000828 -.0029 -.002 -.0015]);
%     
% Acl_FSF = Aaug_FSF - Baug_FSF*Kaug_FSF;
% Bcl_FSF = [Faug_FSF Gaug_FSF];
% Ccl_FSF = Caug_FSF;
% Dcl_FSF = [Daug_FSF zeros(3,1)];
% 
% sys_FSF = ss(Acl_FSF, Bcl_FSF, Ccl_FSF, Dcl_FSF);


r_aug_FSF = repmat([0 .5 0 d],length(ts),1);%reference input augmented with disturbance

[y_FSF,~,x_FSF] = lsim(sys_FSF, r_aug_FSF, ts, x0);
%[y_FSF,~,x_FSF] = lsim(sys_FSF, r_aug_FSF, ts, [x0; zeros(3,1)]); %integral control version

%Check actuators
u = zeros(3,length(ts));
for ii = 1:length(ts)
    u(:,ii) = F*[0;.5;0] - K * x_FSF(ii,:)'; %non-integral version
    %u(:,ii) = -Kaug_FSF * x_FSF(ii,:)'; %integral control version
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
sys = ss(A_CL_aug,B_CL_aug_tot,C_CL_aug,D_CL_aug);
[ys,~,xs] = lsim(sys,us,ts,x_cl_aug0);
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
plot(ax,ts,xs(:,4),'LineWidth',2,'Color','k')
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
sys = ss(A_CL_aug,B_CL_aug_tot,C_CL_aug,D_CL_aug);
[ys,~,xs] = lsim(sys,us,ts,x_cl_aug0);
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
plot(ax,ts,xs(:,4),'LineWidth',2,'Color','k')
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
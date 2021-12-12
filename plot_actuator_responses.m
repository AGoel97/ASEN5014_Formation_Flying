function plot_actuator_responses(ts, rs, xs, F, K, u_max, title)
% Arguments:
%   ts: Time history [sec]
%   rs: Reference history [km]
%   xs: State history [km, km/s]
%   F: Feedforward gain matrix (set to 0 for integral control)
%   K: Feedback gain matrix
%   u_max: Maximum allowable actuator input [km/s^2]
%   name: Plot name

% Reconstruct inputs
u = F*rs - K*xs(:,1:6)';
u = u * 1e3; % Convert to m/s^2
u_max = u_max * 1e3;

figure()
fig = gcf;
fig.Position = [0 50 1000 650];

subplot(3,1,1)
plot(ts,u(1,:),'LineWidth',2,'Color','r')
hold on
plot([0 max(ts)],[u_max u_max],'k:')
plot([0 max(ts)],[-u_max -u_max],'k:')
ylabel('Radial Actuator Response (m/s^2)')

subplot(3,1,2)
plot(ts,u(2,:),'LineWidth',2,'Color','k')
hold on
plot([0 max(ts)],[u_max u_max],'k:')
plot([0 max(ts)],[-u_max -u_max],'k:')
ylabel('In-Track Actuator Response (m/s^2)')

subplot(3,1,3)
plot(ts,u(3,:),'LineWidth',2,'Color','b')
hold on
plot([0 max(ts)],[u_max u_max],'k:')
plot([0 max(ts)],[-u_max -u_max],'k:')
ylabel('Cross-Track Actuator Response (m/s^2)')
xlabel('Time (sec)')

sgtitle(title);

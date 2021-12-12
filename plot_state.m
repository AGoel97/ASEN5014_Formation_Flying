function plot_state(ts, xs, name)
% Arguments:
%   ts: Time history [sec]
%   xs: State history [km, km/s]
%   name: Plot name

figure()
fig = gcf;
fig.Position = [0 50 1000 650];

ax = subplot(3,2,1);
plot(ax,ts,xs(:,1),'LineWidth',2,'Color','r')
hold on;
yline(.25,':');yline(-.25,':');
yline(.05,'--');yline(-.05,'--');
ylabel('x - radial position (km)')
grid on
ax = subplot(3,2,3);
plot(ax,ts,xs(:,2),'LineWidth',2,'Color','k')
hold on;
yline(.5+.25,':');yline(.5-.25,':');
yline(.5+.05,'--');yline(.5-.05,'--');
ylabel('y - along-track position (km)')
grid on
ax = subplot(3,2,5);
plot(ax,ts,xs(:,3),'LineWidth',2,'Color','b')
hold on;
yline(.25,':');yline(-.25,':');
yline(.05,'--');yline(-.05,'--');
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

sgtitle(sprintf('Simulated States (%s)', name));
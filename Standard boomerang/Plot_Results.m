%% Plot and compare aerodynamic forces and moments
close all
load('FA-MA-noapprox.mat')
figure(1)
plot(t,FA(1,:), sim.tout, sim.FA(1,:), 'LineWidth',1.5)
xlabel('Time [s]')
ylabel('Force [N]')
title('F^A_1')
legend('Complete model', 'Approximated model')
xlim([0 min([sim.tout(end) t(end)])])

figure(2)
plot(t,FA(2,:), sim.tout, sim.FA(2,:), 'LineWidth',1.5)
xlabel('Time [s]')
ylabel('Force [N]')
title('F^A_2')
legend('Complete model', 'Approximated model')
xlim([0 min([sim.tout(end) t(end)])])

figure(3)
plot(t,FA(3,:), sim.tout, sim.FA(3,:), 'LineWidth',1.5)
xlabel('Time [s]')
ylabel('Force [N]')
title('F^A_3')
legend('Complete model', 'Approximated model')
xlim([0 min([sim.tout(end) t(end)])])

figure(4)
plot(t,MA(1,:), sim.tout, sim.MA(1,:), 'LineWidth',1.5)
xlabel('Time [s]')
ylabel('Moment [Nm]')
title('M^A_1')
legend('Complete model', 'Approximated model')
xlim([0 min([sim.tout(end) t(end)])])

figure(5)
plot(t,MA(2,:), sim.tout, sim.MA(2,:), 'LineWidth',1.5)
xlabel('Time [s]')
ylabel('Moment [Nm]')
title('M^A_2')
legend('Complete model', 'Approximated model')
xlim([0 min([sim.tout(end) t(end)])])

figure(6)
plot(t,MA(3,:), sim.tout, sim.MA(3,:), 'LineWidth',1.5)
xlabel('Time [s]')
ylabel('Moment [Nm]')
title('M^A_3')
legend('Complete model', 'Approximated model')
xlim([0 min([sim.tout(end) t(end)])])

% figure(7)
% plot(sim.tout, abs(FA(3,1:min([length(t), length(sim.tout)])) - sim.FA(3,1:min([length(t), length(sim.tout)])) )./FA(3,1:min([length(t), length(sim.tout)])) )

mean(abs(FA(3,1:min([length(t), length(sim.tout)])) - sim.FA(3,1:min([length(t), length(sim.tout)])) )./FA(3,1:min([length(t), length(sim.tout)])))

mean(abs(MA(1,1:min([length(t), length(sim.tout)])) - sim.MA(1,1:min([length(t), length(sim.tout)])) )./abs((MA(1,1:min([length(t), length(sim.tout)])))))

%% Plot results over existing plots - Dashed line (comment the 'close all' on main)
R_inertial = zeros(length(sim.tout),3);
u = zeros(length(sim.tout),3);
normU = zeros(length(sim.tout),1);

for i = 1:length(sim.tout)
    R_inertial(i, :) = sim.R_inertial(:,1,i);
    u(i,:) = sim.U(:,1,i);
    normU(i) = norm(u(i,:));
end

% Plots
figure(1);
hold on
plot(sim.tout, R_inertial(:,3), 'k--')
xlabel('Time [s]')
ylabel('Altitude [m]')

figure(2);
hold on
plot(R_inertial(1,1),R_inertial(1,2), 'go', ...
    R_inertial(:,1), R_inertial(:,2), 'k--', ...
    R_inertial(end,1), R_inertial(end,2), 'ro')
xlabel('Displacement along x [m]')
ylabel('Displacement along y [m]')
title('Trajectory - top view (xy plane)')
% legend('', 'Pitching moment$\ne$0', '','','Pitching moment=0', '','Interpreter','latex')

figure(3);
hold on
plot(R_inertial(1,1),R_inertial(1,3), 'go', ...
    R_inertial(:,1), R_inertial(:,3), 'k--', ...
    R_inertial(end,1), R_inertial(end,3), 'ro')
xlabel('Displacement along x [m]')
ylabel('Displacement along z [m]')
title('Trajectory - xz plane')
% legend('', 'Pitching moment$\ne$0', '','','Pitching moment=0', '','Interpreter','latex')

figure(4);
hold on
plot(R_inertial(1,2), R_inertial(1,3), 'go', ...
    R_inertial(:,2), R_inertial(:,3), 'k--', ...
    R_inertial(end,2), R_inertial(end,3), 'ro')
xlabel('Displacement along y [m]')
ylabel('Displacement along z [m]')
title('Trajectory - yz plane')
% legend('', 'Pitching moment$\ne$0', '','','Pitching moment=0', '','Interpreter','latex')

figure(5);
hold on
plot(sim.tout, normU, 'k--')
xlabel('Time [s]')
ylabel('Speed U [m/s]')

figure(6);
hold on
plot3(R_inertial(1,1), R_inertial(1,2), R_inertial(1,3), 'go', ...
    R_inertial(:,1), R_inertial(:,2), R_inertial(:,3), 'k--', ...
    R_inertial(end,1), R_inertial(end,2), R_inertial(end,3), 'ro')
title('3D Trajectory')
%% Plot results over existing plots - Dash-dot line (comment the 'close all' on main)
% Alternative mark for trajectory
R_inertial = zeros(length(sim.tout),3);
u = zeros(length(sim.tout),3);
normU = zeros(length(sim.tout),1);

for i = 1:length(sim.tout)
    R_inertial(i, :) = sim.R_inertial(:,1,i);
    u(i,:) = sim.U(:,1,i);
    normU(i) = norm(u(i,:));
end

% Plots
figure(1);
hold on
plot(sim.tout, R_inertial(:,3), 'k-.')
xlabel('Time [s]')
ylabel('Altitude [m]')
% legend('\phi_0 = 70°', '\phi_0 = 80°', '\phi_0 = 60°')

figure(2);
hold on
plot(R_inertial(1,1),R_inertial(1,2), 'go', ...
    R_inertial(:,1), R_inertial(:,2), 'k-.', ...
    R_inertial(end,1), R_inertial(end,2), 'ro')
xlabel('Displacement along x [m]')
ylabel('Displacement along y [m]')
title('Trajectory - top view (xy plane)')
% legend('Starting point', '\phi_0 = 70°', '', '', '\phi_0 = 80°', '', '', '\phi_0 = 60°', 'End point')

figure(3);
hold on
plot(R_inertial(1,1),R_inertial(1,3), 'go', ...
    R_inertial(:,1), R_inertial(:,3), 'k-.', ...
    R_inertial(end,1), R_inertial(end,3), 'ro')
xlabel('Displacement along x [m]')
ylabel('Displacement along z [m]')
title('Trajectory - xz plane')
% legend('Starting point', '\phi_0 = 70°', '', '', '\phi_0 = 80°', '', '', '\phi_0 = 60°', 'End point')

figure(4);
hold on
plot(R_inertial(1,2), R_inertial(1,3), 'go', ...
    R_inertial(:,2), R_inertial(:,3), 'k-.', ...
    R_inertial(end,2), R_inertial(end,3), 'ro')
xlabel('Displacement along y [m]')
ylabel('Displacement along z [m]')
title('Trajectory - yz plane')
% legend('Starting point', '\phi_0 = 70°', '', '', '\phi_0 = 80°', '', '', '\phi_0 = 60°', 'End point')

figure(5);
hold on
plot(sim.tout, normU, 'k-.')
xlabel('Time [s]')
ylabel('Speed U [m/s]')
% legend('\phi_0 = 70°', '\phi_0 = 80°', '\phi_0 = 60°')

figure(6);
hold on
plot3(R_inertial(1,1), R_inertial(1,2), R_inertial(1,3), 'go', ...
    R_inertial(:,1), R_inertial(:,2), R_inertial(:,3), 'k-.', ...
    R_inertial(end,1), R_inertial(end,2), R_inertial(end,3), 'ro')
title('3D Trajectory')
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
% legend('Starting point', '\phi_0 = 70°', '', '', '\phi_0 = 80°', '', '', '\phi_0 = 60°', 'End point')

%% Plot the reference trajectory
load('Returning_Traj_XYZ-U-tout')

% Plots
figure(1)
plot(t_sim1, R_inertial1(:,3), 'k-')
xlabel('Time [s]')
ylabel('Altitude [m]')

figure(2)
plot(R_inertial1(1,1),R_inertial1(1,2), 'go', ...
    R_inertial1(:,1), R_inertial1(:,2), 'k-', ...
    R_inertial1(end,1), R_inertial1(end,2), 'ro')
xlabel('Displacement along x [m]')
ylabel('Displacement along y [m]')
title('Trajectory - top view (xy plane)')
axis equal
% xlim([-10 10])
% ylim([-5 18])

figure(3)
plot(R_inertial1(1,1),R_inertial1(1,3), 'go', ...
    R_inertial1(:,1), R_inertial1(:,3), 'k-', ...
    R_inertial1(end,1), R_inertial1(end,3), 'ro')
xlabel('Displacement along x [m]')
ylabel('Displacement along z [m]')
title('Trajectory - xz plane')
axis equal
% xlim([-10 10])
% ylim([-2 8])

figure(4)
plot(R_inertial1(1,2), R_inertial1(1,3), 'go', ...
    R_inertial1(:,2), R_inertial1(:,3), 'k-', ...
    R_inertial1(end,2), R_inertial1(end,3), 'ro')
xlabel('Displacement along y [m]')
ylabel('Displacement along z [m]')
title('Trajectory - yz plane')
axis equal
% xlim([-5 18])
% ylim([-2 8])

figure(5)
plot(t_sim1, normU1, 'k-')
xlabel('Time [s]')
ylabel('Speed U [m/s]')
% xlim([0 5])
% ylim([0 25])

figure(6)
plot3(R_inertial1(1,1), R_inertial1(1,2), R_inertial1(1,3), 'go', ...
    R_inertial1(:,1), R_inertial1(:,2), R_inertial1(:,3), 'k-', ...
    R_inertial1(end,1), R_inertial1(end,2), R_inertial1(end,3), 'ro')
title('3D Trajectory')
grid on

%% Comparison with Vassberg model
close all
Cl0 = 0.15;         % Approx coefficient of lift at 0 deg AoA
Rp = 3*m*sin(PHI0)/(rho*S*Cl0);

f = @(x,y) (x-Rp*cosd(ThAng + 90)).^2 + (y - Rp*sind(ThAng + 90)).^2 - Rp^2;

fimplicit(f, [-2*Rp 2*Rp])
figure(1);
hold on
plot(R_inertial(1,1),R_inertial(1,2), 'go', ...
    R_inertial(:,1), R_inertial(:,2), 'k-', ...
    R_inertial(end,1), R_inertial(end,2), 'ro')
xlabel('Displacement along x [m]')
ylabel('Displacement along y [m]')
legend('Vassberg model', '', 'Computed trajectory', '')
axis equal

%%
Cl0 = 0.15;         % Approx coefficient of lift at 0 deg AoA
% Starting point
x1 = R_pos0(1);
y1 = R_pos0(2);

% Desired point to reach
x2 = 5;
y2 = 2;

% Approximate boomerang flight path radius (works well for half a
% revolution)

Rp = 3*m*sin(PHI0)/(rho*S*Cl0);

% Finding possible centers of the circumference
syms a b
[a,b] = solve((a-x1)^2+(b-y1)^2 == Rp^2,(a-x2)^2+(b-y2)^2==Rp^2,a,b);

f = @(x) -sqrt(abs(Rp^2 - (x-a(2)).^2)) + b(2);

xplot = linspace(min(x1,x2), max(x1,x2), 1000);
yplot = f(xplot);

plot(xplot,yplot)

%plot arc
% syms x y
% fimplicit((x-a(1))^2 + (y-b(1))^2 == Rp^2, [min(x1,x2),max(x1,x2), ...
%     min(y1,y2),max(y1,y2)])
% axis equal
% 
% figure
% fimplicit((x-a(2))^2 + (y-b(2))^2 == Rp^2, [min(x1,x2),max(x1,x2), ...
%     min(y1,y2),max(y1,y2)])
% figure
% fimplicit((x-a(2))^2 + (y-b(2))^2 == Rp^2, [-2*Rp 2*Rp -2*Rp 2*Rp])
% axis equal
%%
%first input
a=[0 1]; %P1
b=[1 0]; %P2
r=1;     %radius
%next solution
syms x y
[x,y]=solve((x-a(1))^2+(y-a(2))^2==r^2,(x-b(1))^2+(y-b(2))^2==r^2,x,y);
%plot arc
syms x y
ezplot((x-x(1))^2+(y-y(1))^2==r^2,[min(a(1),b(1)),max(a(1),b(1)), ...
    min(a(2),b(2)),max(a(2),b(2))])
axis equal
figure
ezplot((x-x(2))^2+(y-y(2))^2==r^2,[min(a(1),b(1)),max(a(1),b(1)), ...
    min(a(2),b(2)),max(a(2),b(2))])
axis equal

%% Plots 
load('Vassberg.mat')
ux = s(:,1);
uy = s(:,2);
uz = s(:,3);
u = [ux, uy, uz];

phi = s(:,4);
theta = s(:,5);
psi = s(:,6);
Rp = 3*m*sin(PHI0)/(rho*S*Cl0);  % Flight path radius
% circle = @(x,y) x^2 + (y-Rp)^2 - Rp^2;
POS = zeros(length(t), 3);
POS(1,:) = R_pos0;

U = zeros(length(t),3);

for i = 1:length(t)
    R01 = [
    1, 0, 0;
    0, cos(phi(i)), sin(phi(i));
    0, -sin(phi(i)), cos(phi(i))
    ];
    R02 = [
    cos(theta(i)), 0, -sin(theta(i));
    0, 1, 0;
    sin(theta(i)), 0, cos(theta(i))
    ];
    R03 = [
    cos(psi(i)), sin(psi(i)), 0;
    -sin(psi(i)), cos(psi(i)), 0;
    0, 0, 1
    ];
    T0 = R01*R02*R03;
    invT0 = transpose(T0);
    U(i,:) = invT0*u(i,:)';
end

POSx = cumtrapz(t,U(:,1)) + POS(1,1);
POSy = cumtrapz(t,U(:,2)) + POS(1,2);
POSz = cumtrapz(t,U(:,3)) + POS(1,3);

cut = POSz > 0;

POS = [POSx(cut) POSy(cut) POSz(cut)];


figure(1);
hold on
plot(t(cut),POS(:,3), '--k')
title('Altitude over time')
xlabel('Time [s]')
ylabel('Altitude [s]')

figure(2);
hold on
plot(POS(1,1), POS(1,2), 'go', ...
    POS(:,1), POS(:,2), 'k--', ...
    POS(end,1), POS(end,2), 'ro')
xlabel('Displacement along x [m]')
ylabel('Displacement along y [m]')
title('Trajectory - top view (xy plane)')
legend('', 'Reference model', '','','Vassberg model', '','Interpreter','latex')
% hold on
% fimplicit(circle)
% axis equal

figure(3);
hold on
plot(POS(1,1),POS(1,3), 'go', ...
    POS(:,1), POS(:,3), 'k--', ...
    POS(end,1), POS(end,3), 'ro')
xlabel('Displacement along x [m]')
ylabel('Displacement along z [m]')
title('Trajectory - xz plane')
legend('', 'Reference model', '','','Vassberg model', '','Interpreter','latex')
% axis equal

figure(4);
hold on
plot(POS(1,2), POS(1,3), 'go', ...
    POS(:,2), POS(:,3), 'k--', ...
    POS(end,2), POS(end,3), 'ro')
xlabel('Displacement along y [m]')
ylabel('Displacement along z [m]')
title('Trajectory - yz plane')
legend('', 'Reference model', '','','Vassberg model', '','Interpreter','latex')
% axis equal

figure(5);
hold on
plot(t(cut), vecnorm(u(cut,:)'))

% figure(6)
% plot3(POS(1,1), POS(1,2), POS(1,3), 'go', ...
%     POS(:,1), POS(:,2), POS(:,3),'k--', ...
%     POS(end,1), POS(end,2), POS(end,3), 'ro')
% grid on
% xlabel('X axis')
% ylabel('Y axis')
% zlabel('Z axis')

%% Plot top view trajectory only - PID vs PD
R_inertial = zeros(length(sim.tout),3);
u = zeros(length(sim.tout),3);
normU = zeros(length(sim.tout),1);

%-Load the correct file for uncontrolled trajectory under effects of wind-%
load('Wind_wx_Traj_XYZ-U-tout.mat')
load('Pitch_control_trajectory_2wx.mat')

x_w = R_inertial_w(:,1);
y_w = R_inertial_w(:,2);

z_und = R_inertial1(:,3);

R_inertial(i, :) = sim.R_inertial(:,:);
% u(i,:) = sim.U(:,:);
% normU = vecnorm(u);


figure(1)
plot(R_inertial(1,1),R_inertial(1,2), 'go', ...
    R_inertial(:,1), R_inertial(:,2), 'k-', ...
    R_inertial(end,1), R_inertial(end,2), 'ro')
xlabel('Displacement along x [m]')
ylabel('Displacement along y [m]')
title('Trajectory - top view (xy plane)')
hold on
plot(R_inertial_c(:,1), R_inertial_c(:,2), 'k-.')
hold on
plot(x_des,y_des, 'b-');
hold on
plot(x_w, y_w, 'r--');
legend('', 'PID trajectory', '', 'PD trajectory', 'Desired trajectory','Uncontrolled trajectory', 'Location','southeast')
axis equal

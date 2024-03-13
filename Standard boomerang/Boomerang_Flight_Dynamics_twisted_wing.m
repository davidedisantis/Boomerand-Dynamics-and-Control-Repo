%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                   Boomerang trajectory simulation                   %%%
%%%                   Davide Di Santis - November 2023                  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Source: Flight Dynamics of the Boomerang, Part 1: Fundamental Analysis%%
%% Data
clear
% close all
clc

tic
% Simulation parameters
t = 6;                     % Simulation time [s]
dt = 1e-3;                 % Fixed-step timestep
minstep = 1e-5;            % Min timestep for variable-step simulation [s]
maxstep = 1.66e-4;         % Max timestep for variable-step simulation [s]

t_action1 = 1;             % Instant at which a change of wing pitch occurs
t_action2 = 2;             % Instant at which a second change of pitch occurs
delta_theta1 = 0;          % Pitch variation at instant t_action1 [deg]
delta_theta2 = 0;          % Pitch variation at instant t_action2 (relative to last pitch variation) [deg]
t_rise = 0.01;             % Time needed to rotate by 1 deg [s/deg]
t_rise = t_rise*180/pi;    % [s/rad]
delta_theta1 = deg2rad(delta_theta1);
delta_theta2 = deg2rad(delta_theta2);

% Planetary parameters
g = 9.81;                % Gravity acceleration [m/s^2]
rho = 1.22;              % Air density at sea level [kg/m^3]
v_wind = [0; 0; 0];      % Wind speed in inertial reference frame [m/s]

% Boomerang design parameters
l_blade = 3e-1; % Length of one blade [m]
nw = 2;         % Number of wings [ ]
R = 2.69e-1;    % Radius of rotation [m]
S = 2.28e-1;    % Disk area [m^2]
m = 1.30e-1;    % Boomerang mass [kg]
c = 4.88e-2;    % Mean chord [m]

x_ac = 6.1e-2; % Position of aerodynamic center in body coordinates
LAMBDAj = 120;  % Wing sweep angle [deg]
LAMBDAj = deg2rad(LAMBDAj);
gamma = 120;    % Folding angle [deg]
gamma = deg2rad(gamma);
betaj = 0;      % Wing coning angle 
betaj = deg2rad(betaj);
thetaj0 = 3;     % Wing pitch angle at root
thetaj0 = deg2rad(thetaj0);
k = -0.1745;    % thetaj = thetaj0 + k*wingspan

% Moments of inertia of a single blade
I_xi = 1.88e-3; 
I_eta = 4.88e-6;
I_zeta = 1.92e-3;
I_xieta = 5.14e-21;

Jj = [
    I_xi I_xieta 0;
    I_xieta I_eta 0;
    0 0 I_zeta
    ];

% Launch conditions
ThAng = 45;                                 % Throw angle from East direction
U0 = [22*cosd(ThAng); 22*sind(ThAng); 10];   % Initial throw speed in inertial frame [m/s]
omega0 = [0; 0; 10*pi*2];                   % Initial angular speed [rad/s]
PHI0 = 77;                                  % Initial roll angle (Between non-spinning frame and inertial) [deg]
PHI0 = deg2rad(PHI0);
PSI0 = 240;                                 % Initial yaw angle (Between non-spinning frame and inertial) [deg]
PSI0 = deg2rad(PSI0);
THETA0 = 0;                                 % Initial pitch angle (Between non-spinning frame and inertial) [deg]
THETA0 = deg2rad(THETA0);  
R_pos0 = [0; 0; 1.5];                       % Initial position in inertial frame

% Parameters for integration of aerodynamic forces
n = 50; % Number of intervals
l_integrate = linspace(0,l_blade,n+1);
h = l_blade/n;

% Rotation matrices
Tj = zeros(3,3, nw);
invTj = zeros(3,3,nw);

for i = 0:(nw-1)
    Rj1 = [
    1, 0, 0;
    0, cos(betaj), sin(betaj);
    0, -sin(betaj), cos(betaj)
    ];
    Rj2 = [
    cos(thetaj0), 0, -sin(thetaj0);
    0, 1, 0;
    sin(thetaj0), 0, cos(thetaj0)
    ];
    Rj3 = [
    cos(LAMBDAj-pi/2 + gamma*i), sin(LAMBDAj-pi/2 + gamma*i), 0;
    -sin(LAMBDAj-pi/2 + gamma*i), cos(LAMBDAj-pi/2+ gamma*i), 0;
    0, 0, 1
    ];
    Tj(:,:,i+1) = Rj2*Rj1*Rj3;
    % Tj(2,2,i+1) = -Tj(2,2,i+1);
    invTj(:,:,i+1) = transpose(Tj(:,:,i+1));
end

% Computation of moment of inertia wrt body frame
Ji = zeros(3,3,nw);
J = zeros(3);

for i = 1:nw
    Ji(:,:,i) = invTj(:,:,i)*Jj*Tj(:,:,i);
    J = J + (Ji(:,:,i));
end
J1 = invTj(:,:,1)*Jj*Tj(:,:,1);
J = J1*2;

J(2,2) = J(2,2) + m*x_ac^2;
J(3,3) = J(3,3) + m*x_ac^2;

J = diag(diag(J));

invJ = inv(J);

RI01 = [
    1, 0, 0;
    0, cos(PHI0), sin(PHI0);
    0, -sin(PHI0), cos(PHI0)
    ];
RI02 = [
    cos(THETA0), 0, -sin(THETA0);
    0, 1, 0;
    sin(THETA0), 0, cos(THETA0)
];
RI03 = [
    cos(PSI0), sin(PSI0), 0;
    -sin(PSI0), cos(PSI0), 0;
    0, 0, 1
];
TI0 = RI01*RI02*RI03;
invTI0 = transpose(TI0);

% Starting conditions, body frame
u0 = TI0*U0;
r0 = TI0*R_pos0;

Eul0 = [0; 0; 0];
% Quat0 = [1; 0; 0; 0];

% CL data
M1 = readmatrix("Cl.csv");
x_CL = M1(:,3);
y_CL = 1.23*M1(:,2);
% y_CL = M1(:,2);
CLdata = [x_CL, y_CL];

% CD data
M2 = readmatrix("Cd.csv");
x_CD = M2(:,3);
y_CD = 1.03*M2(:,2);
% y_CD = M2(:,2);
CDdata = [x_CD, y_CD];

% CM data
M3 = readmatrix("Cm.csv");
x_CM = M3(:,3);
y_CM = M3(:,2);
CMdata = [x_CM, y_CM];

sim = sim("Boomerang_twisted_wing_Simulink.slx");

toc
%% Plot results 

phi = sim.EulAng(:,1);
theta = sim.EulAng(:,2);
psi = sim.EulAng(:,3);
R_inertial = sim.R_inertial(:,:)';
u = sim.U(:,:)';
normU = vecnorm(u');

% Plots
figure(1)
plot(sim.tout, R_inertial(:,3), 'k-')
xlabel('Time [s]')
ylabel('Altitude [m]')

figure(2)
plot(R_inertial(1,1),R_inertial(1,2), 'go', ...
    R_inertial(:,1), R_inertial(:,2), 'k-', ...
    R_inertial(end,1), R_inertial(end,2), 'ro')
xlabel('Displacement along x [m]')
ylabel('Displacement along y [m]')
title('Trajectory - top view (xy plane)')
axis equal
% xlim([-10 10])
% ylim([-5 18])

figure(3)
plot(R_inertial(1,1),R_inertial(1,3), 'go', ...
    R_inertial(:,1), R_inertial(:,3), 'k-', ...
    R_inertial(end,1), R_inertial(end,3), 'ro')
xlabel('Displacement along x [m]')
ylabel('Displacement along z [m]')
title('Trajectory - xz plane')
axis equal
% xlim([-10 10])
% ylim([-2 8])

figure(4)
plot(R_inertial(1,2), R_inertial(1,3), 'go', ...
    R_inertial(:,2), R_inertial(:,3), 'k-', ...
    R_inertial(end,2), R_inertial(end,3), 'ro')
xlabel('Displacement along y [m]')
ylabel('Displacement along z [m]')
title('Trajectory - yz plane')
axis equal
% xlim([-5 18])
% ylim([-2 8])

figure(5)
plot(sim.tout, normU, 'k-')
xlabel('Time [s]')
ylabel('Speed U [m/s]')
% xlim([0 5])
% ylim([0 25])

figure(6)
plot3(R_inertial(1,1), R_inertial(1,2), R_inertial(1,3), 'go', ...
    R_inertial(:,1), R_inertial(:,2), R_inertial(:,3), 'k-', ...
    R_inertial(end,1), R_inertial(end,2), R_inertial(end,3), 'ro')
title('3D Trajectory')
grid on
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

figure(3);
hold on
plot(R_inertial(1,1),R_inertial(1,3), 'go', ...
    R_inertial(:,1), R_inertial(:,3), 'k--', ...
    R_inertial(end,1), R_inertial(end,3), 'ro')
xlabel('Displacement along x [m]')
ylabel('Displacement along z [m]')
title('Trajectory - xz plane')

figure(4);
hold on
plot(R_inertial(1,2), R_inertial(1,3), 'go', ...
    R_inertial(:,2), R_inertial(:,3), 'k--', ...
    R_inertial(end,2), R_inertial(end,3), 'ro')
xlabel('Displacement along y [m]')
ylabel('Displacement along z [m]')
title('Trajectory - yz plane')

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
plot(t_vec, R_inertial1(:,3), 'k-')
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
xlim([-10 10])
ylim([-5 18])

figure(3)
plot(R_inertial1(1,1),R_inertial1(1,3), 'go', ...
    R_inertial1(:,1), R_inertial1(:,3), 'k-', ...
    R_inertial1(end,1), R_inertial1(end,3), 'ro')
xlabel('Displacement along x [m]')
ylabel('Displacement along z [m]')
title('Trajectory - xz plane')
axis equal
xlim([-10 10])
ylim([-2 8])

figure(4)
plot(R_inertial1(1,2), R_inertial1(1,3), 'go', ...
    R_inertial1(:,2), R_inertial1(:,3), 'k-', ...
    R_inertial1(end,2), R_inertial1(end,3), 'ro')
xlabel('Displacement along y [m]')
ylabel('Displacement along z [m]')
title('Trajectory - yz plane')
axis equal
xlim([-5 18])
ylim([-2 8])

figure(5)
plot(t_vec, normU1, 'k-')
xlabel('Time [s]')
ylabel('Speed U [m/s]')
xlim([0 5])
ylim([0 25])

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
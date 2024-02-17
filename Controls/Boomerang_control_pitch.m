%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                    Boomerang trajectory control                     %%%
%%%                   Davide Di Santis - January 2024                   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Source: Flight Dynamics of the Boomerang, Part 1: Fundamental Analysis%%
%% Data
clear
% close all
clc

tic
% Simulation parameters
t = 3.8;                     % Simulation time [s]
dt = 1e-3;                 % Timestep for fixed-step simulation [s]
minstep = 1e-5;            % Min timestep for variable-step simulation [s]
maxstep = 1.66e-4;         % Max timestep for variable-step simulation [s]
servo_rpm = 500;           % Rotation speed of servo
servo_rpm = deg2rad(servo_rpm);

% Desired trajectory: returning path 
load('Returning_Traj_XYZ-U-tout.mat');
x_des = R_inertial1(:,1);
y_des = R_inertial1(:,2);

X_des = [x_des, y_des];

% % Desired trajectory: circular
% R_c = 5;                                % Radius of desired circular trajectory [m]
% 
% x_des1 = 0:dt*4:R_c;
% y_des1 = -sqrt((R_c^2 - x_des1.^2)) + R_c;
% 
% x_des2 = 0:dt*4:R_c;
% y_des2 = sqrt((R_c^2 - x_des2.^2)) + R_c;
% 
% x_des3 = -R_c:dt*4:0;
% y_des3 = sqrt((R_c^2 - x_des3.^2)) + R_c;
% 
% x_des4 = -R_c:dt*4:0;
% y_des4 = -sqrt((R_c^2 - x_des4.^2)) + R_c;
% 
% X_des = [x_des1, x_des2, x_des3, x_des4; y_des1, y_des2, y_des3, y_des4]';

% PID parameters
Kp = 1; 
Ki = 1.2;
Kd = 0.8;

% Planetary parameters
g = 9.81;               % Gravity acceleration [m/s^2]
rho = 1.22;             % Air density at sea level [kg/m^3]
v_wind = [-3; 0; 0];  % Wind velocity [m/s]

% Boomerang design parameters
l_blade = 3e-1; % Length of one blade [m]
nw = 2;         % Number of wings [ ]
R = 2.69e-1;    % Radius of rotation [m]
S = 2.28e-1;    % Disk area [m^2]
m = 1.30e-1;    % Boomerang mass [kg]
c = 4.88e-2;    % Mean chord [m]

x_ac = 7.45e-2; % Position of aerodynamic center in body coordinates
LAMBDAj = 120;  % Wing sweep angle [deg]
LAMBDAj = deg2rad(LAMBDAj);
gamma = 120;    % Folding angle [deg]
gamma = deg2rad(gamma);
betaj = 0;      % Wing coning angle 
betaj = deg2rad(betaj);
thetaj = 0;     % Wing pitch angle
thetaj = deg2rad(thetaj);

% Moments of inertia of a single blade
I_xi = 1.88e-3; 
I_eta = 4.78e-6;
I_zeta = 1.95e-3;
I_xieta = 5.14e-21;

Jj = [
    I_xi I_xieta 0;
    I_xieta I_eta 0;
    0 0 I_zeta
    ];

% Launch conditions
U0 = [25*cosd(30); 25*sind(30); 0];  % Initial throw speed in inertial frame [m/s]
omega0 = [0; 0; 10*pi*2];            % Initial angular speed [rad/s]
PHI0 = 70;                           % Initial roll angle (Between non-spinning frame and inertial) [deg]
PHI0 = deg2rad(PHI0);
PSI0 = 240;                           % Initial yaw angle (Between non-spinning frame and inertial) [deg]
PSI0 = deg2rad(PSI0);
THETA0 = 0;                          % Initial pitch angle (Between non-spinning frame and inertial) [deg]
THETA0 = deg2rad(THETA0);  
R_pos0 = [0; 0; 5];                % Initial position in inertial frame

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
    cos(thetaj), 0, -sin(thetaj);
    0, 1, 0;
    sin(thetaj), 0, cos(thetaj)
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
y_CL = .82*M1(:,2);
% y_CL = M1(:,2);
CLdata = [x_CL, y_CL];

% CD data
M2 = readmatrix("Cd.csv");
x_CD = M2(:,3);
y_CD = .85*M2(:,2);
% y_CD = M2(:,2);
CDdata = [x_CD, y_CD];

% CM data
M3 = readmatrix("Cm.csv");
x_CM = M3(:,3);
y_CM = M3(:,2);
CMdata = [x_CM, y_CM];

% sim = sim("Azuma_Dynamics_Simulink_EulAng_Jbody.slx");
sim = sim("Boomerang_pitch_simulink.slx");

toc
%% Plot top view trajectory only 
R_inertial = zeros(length(sim.tout),3);
u = zeros(length(sim.tout),3);
normU = zeros(length(sim.tout),1);
load('Wind_2wx_Traj_XYZ-U-tout.mat')
load('Pitch_control_trajectory_2wx.mat')

x_w = R_inertial_wx(:,1);
y_w = R_inertial_wx(:,2);

z_und = R_inertial1(:,3) + 3.5;
for i = 1:length(sim.tout)
    R_inertial(i, :) = sim.R_inertial(:,1,i);
    u(i,:) = sim.U(:,1,i);
    normU(i) = norm(u(i,:));
end

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

% figure(1)
% plot(R_inertial(1,1),R_inertial(1,3), 'go', ...
%     R_inertial(:,1), R_inertial(:,3), 'k-', ...
%     R_inertial(end,1), R_inertial(end,3), 'ro')
% xlabel('Displacement along x [m]')
% ylabel('Displacement along z [m]')
% title('Trajectory - back view (xz plane)')
% hold on
% plot(x_des, z_und, 'b-');
% legend('', 'Actual trajectory', '', 'Undisturbed trajectory', 'Location','southeast')
% axis equal

% figure(2)
% plot(R_inertial(1,2),R_inertial(1,3), 'go', ...
%     R_inertial(:,2), R_inertial(:,3), 'k-', ...
%     R_inertial(end,2), R_inertial(end,3), 'ro')
% xlabel('Displacement along y [m]')
% ylabel('Displacement along z [m]')
% title('Trajectory - side view (yz plane)')
% hold on
% plot(y_des, z_und, 'b-');
% legend('', 'Actual trajectory', '', 'Undisturbed trajectory', 'Location','southeast')
% axis equal
% 
% figure(3)
% plot(R_inertial(1,1),R_inertial(1,2), 'go', ...
%     R_inertial(:,1), R_inertial(:,2), 'k-', ...
%     R_inertial(end,1), R_inertial(end,2), 'ro')
% xlabel('Displacement along x [m]')
% ylabel('Displacement along y [m]')
% title('Trajectory - top view (xy plane)')
% hold on
% plot(x_des,y_des, 'b');
% hold on
% plot(x_w, y_w, 'r--');
% legend('', 'Actual trajectory', '', 'Desired trajectory', 'Uncontrolled trajectory', 'Location','southeast')
% axis equal
%% Plot and compare energies
load('Uncontrolled_Energies.mat')
figure(1)
plot(sim.P_E, 'k--')
hold on
plot(P_Enc + 3.5*m*g, 'k') % Adding 3.5*m*g to shift the curve as starting at 5 m 
title('Potential Energy over time')
xlabel('Time [s]')
ylabel('Energy [J]')
legend('E_P no control', 'E_P control')

figure(2)
plot(sim.K_E, 'k--')
hold on
plot(K_Enc, 'k')
title('Kinetic Energy over time')
xlabel('Time [s]')
ylabel('Energy [J]')
legend('E_K no control', 'E_K control')

figure(3)
plot(sim.R_E, 'k--')
hold on
plot(R_Enc, 'k')
title('Rotational Energy over time')
xlabel('Time [s]')
ylabel('Energy [J]')
legend('E_R no control', 'E_R control')
%% Plot top view trajectory only - Desired circular trajectory
R_inertial = zeros(length(sim.tout),3);
u = zeros(length(sim.tout),3);
normU = zeros(length(sim.tout),1);

for i = 1:length(sim.tout)
    R_inertial(i, :) = sim.R_inertial(:,1,i);
    u(i,:) = sim.U(:,1,i);
    normU(i) = norm(u(i,:));
end

figure(1)
plot(R_inertial(1,1),R_inertial(1,2), 'go', ...
    R_inertial(:,1), R_inertial(:,2), 'k-', ...
    R_inertial(end,1), R_inertial(end,2), 'ro')
xlabel('Displacement along x [m]')
ylabel('Displacement along y [m]')
title('Trajectory - top view (xy plane)')
hold on
plot(x_des1,y_des1, 'b')
hold on
plot(x_des2,y_des2, 'b')
hold on
plot(x_des3,y_des3, 'b')
hold on
plot(x_des4,y_des4, 'b')
hold on
plot(R_inertial1(:,1), R_inertial1(:,2), 'r--')
legend('', 'Actual trajectory', '', 'Desired trajectory', '','','', 'Reference trajectory','Location','southeast')
axis equal
%% Plot results 
R_inertial = zeros(length(sim.tout),3);
u = zeros(length(sim.tout),3);
normU = zeros(length(sim.tout),1);

for i = 1:length(sim.tout)
    R_inertial(i, :) = sim.R_inertial(:,1,i);
    u(i,:) = sim.U(:,1,i);
    normU(i) = norm(u(i,:));
end

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
xlim([-10 10])
ylim([-5 18])

figure(3)
plot(R_inertial(1,1),R_inertial(1,3), 'go', ...
    R_inertial(:,1), R_inertial(:,3), 'k-', ...
    R_inertial(end,1), R_inertial(end,3), 'ro')
xlabel('Displacement along x [m]')
ylabel('Displacement along z [m]')
title('Trajectory - xz plane')
axis equal
xlim([-10 10])
ylim([-2 8])

figure(4)
plot(R_inertial(1,2), R_inertial(1,3), 'go', ...
    R_inertial(:,2), R_inertial(:,3), 'k-', ...
    R_inertial(end,2), R_inertial(end,3), 'ro')
xlabel('Displacement along y [m]')
ylabel('Displacement along z [m]')
title('Trajectory - yz plane')
axis equal
xlim([-5 18])
ylim([-2 8])

figure(5)
plot(sim.tout, normU, 'k-')
xlabel('Time [s]')
ylabel('Speed U [m/s]')
xlim([0 5])
ylim([0 25])

figure(6)
plot3(R_inertial(1,1), R_inertial(1,2), R_inertial(1,3), 'go', ...
    R_inertial(:,1), R_inertial(:,2), R_inertial(:,3), 'k-', ...
    R_inertial(end,1), R_inertial(end,2), R_inertial(end,3), 'ro')
title('3D Trajectory')
grid on
%% Plot results over existing plots (comment the 'close all' on main)
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

figure(2);
hold on
plot(R_inertial(1,1),R_inertial(1,2), 'go', ...
    R_inertial(:,1), R_inertial(:,2), 'k-.', ...
    R_inertial(end,1), R_inertial(end,2), 'ro')
xlabel('Displacement along x [m]')
ylabel('Displacement along y [m]')
title('Trajectory - top view (xy plane)')

figure(3);
hold on
plot(R_inertial(1,1),R_inertial(1,3), 'go', ...
    R_inertial(:,1), R_inertial(:,3), 'k-.', ...
    R_inertial(end,1), R_inertial(end,3), 'ro')
xlabel('Displacement along x [m]')
ylabel('Displacement along z [m]')
title('Trajectory - xz plane')

figure(4);
hold on
plot(R_inertial(1,2), R_inertial(1,3), 'go', ...
    R_inertial(:,2), R_inertial(:,3), 'k-.', ...
    R_inertial(end,2), R_inertial(end,3), 'ro')
xlabel('Displacement along y [m]')
ylabel('Displacement along z [m]')
title('Trajectory - yz plane')

figure(5);
hold on
plot(sim.tout, normU, 'k-.')
xlabel('Time [s]')
ylabel('Speed U [m/s]')

figure(6);
hold on
plot3(R_inertial(1,1), R_inertial(1,2), R_inertial(1,3), 'go', ...
    R_inertial(:,1), R_inertial(:,2), R_inertial(:,3), 'k-.', ...
    R_inertial(end,1), R_inertial(end,2), R_inertial(end,3), 'ro')
title('3D Trajectory')
%% Plot results over existing plots (comment the 'close all' on main)
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
plot(sim.tout, R_inertial(:,3), 'k--')
xlabel('Time [s]')
ylabel('Altitude [m]')
legend('\phi_0 = 70°', '\phi_0 = 80°', '\phi_0 = 60°')

figure(2);
hold on
plot(R_inertial(1,1),R_inertial(1,2), 'go', ...
    R_inertial(:,1), R_inertial(:,2), 'k--', ...
    R_inertial(end,1), R_inertial(end,2), 'ro')
xlabel('Displacement along x [m]')
ylabel('Displacement along y [m]')
title('Trajectory - top view (xy plane)')
legend('Starting point', '\phi_0 = 70°', '', '', '\phi_0 = 80°', '', '', '\phi_0 = 60°', 'End point')

figure(3);
hold on
plot(R_inertial(1,1),R_inertial(1,3), 'go', ...
    R_inertial(:,1), R_inertial(:,3), 'k--', ...
    R_inertial(end,1), R_inertial(end,3), 'ro')
xlabel('Displacement along x [m]')
ylabel('Displacement along z [m]')
title('Trajectory - xz plane')
legend('Starting point', '\phi_0 = 70°', '', '', '\phi_0 = 80°', '', '', '\phi_0 = 60°', 'End point')

figure(4);
hold on
plot(R_inertial(1,2), R_inertial(1,3), 'go', ...
    R_inertial(:,2), R_inertial(:,3), 'k--', ...
    R_inertial(end,2), R_inertial(end,3), 'ro')
xlabel('Displacement along y [m]')
ylabel('Displacement along z [m]')
title('Trajectory - yz plane')
legend('Starting point', '\phi_0 = 70°', '', '', '\phi_0 = 80°', '', '', '\phi_0 = 60°', 'End point')

figure(5);
hold on
plot(sim.tout, normU, 'k--')
xlabel('Time [s]')
ylabel('Speed U [m/s]')
legend('\phi_0 = 70°', '\phi_0 = 80°', '\phi_0 = 60°')

figure(6);
hold on
plot3(R_inertial(1,1), R_inertial(1,2), R_inertial(1,3), 'go', ...
    R_inertial(:,1), R_inertial(:,2), R_inertial(:,3), 'k--', ...
    R_inertial(end,1), R_inertial(end,2), R_inertial(end,3), 'ro')
title('3D Trajectory')
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
legend('Starting point', '\phi_0 = 70°', '', '', '\phi_0 = 80°', '', '', '\phi_0 = 60°', 'End point')
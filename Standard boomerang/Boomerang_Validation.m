%% Data
clear
close all
clc

tic
%--------------------------Simulation parameters--------------------------%
t = 5;                     % Simulation time [s]
dt = 2e-4;                % Fixed-step timestep
minstep = 1e-5;            % Min timestep for variable-step simulation [s]
maxstep = 1.66e-4;         % Max timestep for variable-step simulation [s]

% Planetary parameters
g = 9.81;                % Gravity acceleration [m/s^2]
rho = 1.22;              % Air density at sea level [kg/m^3]
v_wind = [0; 0; 0];      % Wind speed in inertial reference frame [m/s]

% Boomerang design parameters
l_blade = 144e-3; % Length of one blade [m]
nw = 3;         % Number of wings [ ]
R = l_blade;    % Radius of rotation [m]
S = pi*R^2;    % Disk area [m^2]
m = 52e-3;    % Boomerang mass [kg]
c = 3.81e-2;    % Mean chord [m]

x_ac = 0; % Position of aerodynamic center in body coordinates
LAMBDAj = 0;  % Wing sweep angle [deg]
LAMBDAj = deg2rad(LAMBDAj);
gamma = 120;    % Folding angle [deg]
gamma = deg2rad(gamma);
betaj = 0;      % Wing coning angle 
betaj = deg2rad(betaj);
thetaj = 0;     % Wing pitch angle
thetaj = deg2rad(thetaj);

%--------------------------Moments of inertia-----------------------------%
Ixx = 228550e-9;
Ixy = -4078.73e-9;
Ixz = 0.36e-9;
Iyx = -4078.73e-9;
Iyy = 213837.98e-9;
Iyz = 0.62e-9;
Izx = 0.36e-9;
Izy = 0.62e-9;
Izz = 442248.23e-9;

J = [Ixx Ixy Ixz;
    Iyx Iyy Iyz;
    Izx Izy Izz];
invJ = inv(J);

%-------------------------Launch conditions-------------------------------%
ThAng = 45;                                 % Throw angle from East direction
U0 = [23*cosd(ThAng); 23*sind(ThAng); 0];   % Initial throw speed in inertial frame [m/s]
r = 590;                                    % Boomerang's RPM
r = r/60*2*pi;                     
omega0 = [0; 0; r];                         % Initial angular speed [rad/s]
PHI0 = 60;                                  % Initial roll angle (Between non-spinning frame and inertial) [deg]
PHI0 = deg2rad(PHI0);
THETA0 = 0;                                 % Initial pitch angle (Between non-spinning frame and inertial) [deg]
THETA0 = deg2rad(THETA0);  
PSI0 = 225;                                 % Initial yaw angle (Between non-spinning frame and inertial) [deg]
PSI0 = deg2rad(PSI0);
R_pos0 = [0; 0; 1.5];                       % Initial position in inertial frame

%------------Parameters for integration of aerodynamic forces-------------%
n = 50; % Number of intervals
n = n +1; 
l_integrate = linspace(0,l_blade,n);
%---------------------------Rotation matrices-----------------------------%
Tj = zeros(3,3, nw);
invTj = zeros(3,3,nw);

for i = 0:(nw-1)
    Rj1 = [1, 0, 0; 0, cos(betaj), sin(betaj); 0, -sin(betaj), cos(betaj)];
    Rj2 = [cos(thetaj), 0, -sin(thetaj); 0, 1, 0; sin(thetaj), 0, cos(thetaj)];
    Rj3 = [cos(LAMBDAj-pi/2 + gamma*i), sin(LAMBDAj-pi/2 + gamma*i), 0; ...
        -sin(LAMBDAj-pi/2 + gamma*i), cos(LAMBDAj-pi/2+ gamma*i), 0;0, 0, 1];
    Tj(:,:,i+1) = Rj2*Rj1*Rj3;
    invTj(:,:,i+1) = transpose(Tj(:,:,i+1));
end
RI01 = [1, 0, 0; 0, cos(PHI0), sin(PHI0); 0, -sin(PHI0), cos(PHI0)];
RI02 = [cos(THETA0), 0, -sin(THETA0); 0, 1, 0; sin(THETA0), 0, cos(THETA0)];
RI03 = [cos(PSI0), sin(PSI0), 0; -sin(PSI0), cos(PSI0), 0; 0, 0, 1];
TI0 = RI01*RI02*RI03;
invTI0 = transpose(TI0);

%-------------------Starting conditions in body frame---------------------%
Eul0 = [0; 0; 0];
u0 = TI0*U0;
r0 = TI0*R_pos0;

%------------------------Aerodynamic coefficients-------------------------%
% CL data
M1 = readmatrix("Cl.csv");
x_CL = M1(:,3);
y_CL = M1(:,2);
% y_CL = M1(:,2);
CLdata = [x_CL, y_CL];

% CD data
M2 = readmatrix("Cd.csv");
x_CD = M2(:,3);
y_CD = M2(:,2);
% y_CD = 0.95*M2(:,2);
CDdata = [x_CD, y_CD];

% CM data
M3 = readmatrix("Cm.csv");
x_CM = M3(:,3);
y_CM = M3(:,2);
CMdata = [x_CM, y_CM];

%--------------------------IMU Sensor modeling----------------------------% Present data = ISM330DHCX IMU
% acc = struct('scalefactor', [], 'bias',[], 'lim', [], 'noise', []);
% gyro = struct('scalefactor', [], 'bias', [], 'gbias', [], 'lim',[], 'noise',[]);
% f = 1000;                            % Sampling frequency [Hz]
% acc.sf = [1 0.005 0.005; 0.005 1 0.005; 0.005 0.005 1]; % Scale factor and cross-couplings []
% acc.bias = [1 1 1].*60e-3*g;
% acc.lim = [-16*g  -16*g  -16*g  16*g 16*g 16*g];
% acc.noise = [80e-6*g/sqrt(f) 80e-6*g/sqrt(f) 80e-6*g/sqrt(f)];
% 
% gyro.sf = [1 0.01 0.01; 0.01 1 0.01; 0.01 0.01 1];
% gbias = deg2rad(3);
% gyro.bias = [gbias gbias gbias];
% gyro.gbias = [deg2rad(0.1) deg2rad(0.1) deg2rad(0.1)];
% glim = deg2rad(4000);
% gyro.lim = [-glim  -glim  -glim  glim glim glim];
% gyro.noise = [deg2rad(7e-3)/sqrt(f) deg2rad(7e-3)/sqrt(f) deg2rad(7e-3)/sqrt(f)];
% 
% noise = [acc.noise gyro.noise];

sim = sim("Boomerang_Simulink.slx");

toc

%% Plot results 

phi = sim.EulAng(:,1);
theta = sim.EulAng(:,2);
psi = sim.EulAng(:,3);
R_inertial = sim.R_inertial(:,:)';
U = sim.U(:,:)';
normU = vecnorm(U');

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
xlim([0 5])
ylim([0 25])
% figure(6)
% plot(sim.tout, sim.omega(:,3),'k')

% figure(6)
% plot3(R_inertial(1,1), R_inertial(1,2), R_inertial(1,3), 'go', ...
%     R_inertial(:,1), R_inertial(:,2), R_inertial(:,3), 'k-', ...
%     R_inertial(end,1), R_inertial(end,2), R_inertial(end,3), 'ro')
% title('3D Trajectory')
% grid on
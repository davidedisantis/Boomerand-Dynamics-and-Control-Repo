%% Data
clear
% close all
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
l_blade_0 = 0.5; % Portion of blade that has been cropped
l_blade = 144e-3; % Length of one blade [m]
nw = 3;         % Number of wings [ ]
R = l_blade;    % Radius of rotation [m]
S = pi*R^2;    % Disk area [m^2]
m_b = 35e-3;    % Boomerang mass [kg]
c = 3.81e-2;    % Mean chord [m]
m_add = zeros(3);
% m_add(3,3) = 8/3*rho*R^3;
m = m_b*eye(3) + m_add;
invm = inv(m);

d_ac = c/4; % Distance between CoM and a.c.
LAMBDAj = 0;  % Wing sweep angle [deg]
LAMBDAj = deg2rad(LAMBDAj);
gamma = 120;    % Folding angle [deg]
gamma = deg2rad(gamma);
betaj = 0;      % Wing coning angle 
betaj = deg2rad(betaj);
thetaj = 0;     % Wing pitch angle
thetaj = deg2rad(thetaj);
x_ac = zeros(3,nw);
for i = 1:nw
    Rac = [cos(gamma*i), sin(gamma*i), 0; -sin(gamma*i), cos(gamma*i), 0;0, 0, 1];
    x_ac(:,i) = Rac*[0; -d_ac; 0];
end
    
%--------------------------Moments of inertia-----------------------------%
Ixx = 206217e-9;
Ixy = 0e-9;
Ixz = -3102.96e-9;
Iyx = Ixy;
Iyy = 213838e-9;
Iyz = -1890e-9;
Izx = Ixz;
Izy = Iyz;
Izz = 412199e-9;

J = [Ixx Ixy Ixz;
    Iyx Iyy Iyz;
    Izx Izy Izz];
invJ = inv(J);

%-------------------------Launch conditions-------------------------------%
ThAng = 90;                                 % Throw angle from East direction
U0 = [18*cosd(ThAng); 18*sind(ThAng); 0];   % Initial throw speed in inertial frame [m/s]
r = 730;                                    % Boomerang's RPM
r = r/60*2*pi;                     
omega0 = [0; 0; r];                         % Initial angular speed [rad/s]
PHI0 = 60;                                  % Initial roll angle (Between non-spinning frame and inertial) [deg]
PHI0 = deg2rad(PHI0);
THETA0 = 0;                                 % Initial pitch angle (Between non-spinning frame and inertial) [deg]
THETA0 = deg2rad(THETA0);  
PSI0 = 270;                                 % Initial yaw angle (Between non-spinning frame and inertial) [deg]
PSI0 = deg2rad(PSI0);
R_pos0 = [0; 0; 1.8];                       % Initial position in inertial frame

%------------Parameters for integration of aerodynamic forces-------------%
n = 50; % Number of intervals
n = n +1; 
l_integrate = linspace(l_blade_0,l_blade,n);
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
y_CL = 1.3*M1(:,2);
% y_CL = M1(:,2);
CLdata = [x_CL, y_CL];

% CD data
M2 = readmatrix("Cd.csv");
x_CD = M2(:,3);
y_CD = 1.3*M2(:,2);
% y_CD = 0.95*M2(:,2);
CDdata = [x_CD, y_CD];

% CM data
M3 = readmatrix("Cm.csv");
x_CM = M3(:,3);
y_CM = M3(:,2);
CMdata = [x_CM, y_CM];

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
hold on
plot(sim.tout, R_inertial(:,3), 'k-')
xlabel('Time [s]')
ylabel('Altitude [m]')

figure(2)
hold on
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
hold on
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
hold on
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
hold on
plot(sim.tout, normU, 'k-')
xlabel('Time [s]')
ylabel('Speed U [m/s]')
xlim([0 5])
ylim([0 25])

figure(6)
plot(sim.tout, sim.omega(:,3),'k')

figure(6)
plot3(R_inertial(1,1), R_inertial(1,2), R_inertial(1,3), 'go', ...
    R_inertial(:,1), R_inertial(:,2), R_inertial(:,3), 'k-', ...
    R_inertial(end,1), R_inertial(end,2), R_inertial(end,3), 'ro')
title('3D Trajectory')
grid on
fprintf('\nFinal deviation from the thrower is %f m\n', norm([R_inertial(end,1),R_inertial(end,2)]))
%% Trajectory and attitude
phi = sim.EulAng(:,1);
theta = sim.EulAng(:,2);
psi = sim.EulAng(:,3);
R_inertial = sim.R_inertial(:,:)';

R_inertial1 = R_inertial;
x_v = [1; 0; 0];
y_v = [0; 1; 0];
z_v = [0; 0; 1];
for i = 1:50:length(phi)
    R1 = [
    1, 0, 0;
    0, cos(phi(i)), sin(phi(i));
    0, -sin(phi(i)), cos(phi(i))
    ];
    R2 = [
    cos(theta(i)), 0, -sin(theta(i));
    0, 1, 0;
    sin(theta(i)), 0, cos(theta(i))
    ];
    R3 = [
    cos(psi(i)), sin(psi(i)), 0;
    -sin(psi(i)), cos(psi(i)), 0;
    0, 0, 1
    ];
    T0 = R1*R2*R3;
    invT0 = transpose(T0);
    attx = invTI0*invT0*x_v;
    atty = invTI0*invT0*y_v;
    attz = invTI0*invT0*z_v;
    % figure(1)
    % quiver3(origin(1),origin(2),origin(3), attx(1), attx(2), attx(3), 'b')
    % hold on
    % quiver3(origin(1),origin(2),origin(3), atty(1), atty(2), atty(3), 'b')
    % hold on
    % quiver3(origin(1),origin(2),origin(3), attz(1), attz(2), attz(3), 'r')
    % hold off

    figure(1);
    plot3(R_inertial1(:,1),R_inertial1(:,2),R_inertial1(:,3),'k-')
    hold on
    quiver3(R_inertial1(i,1),R_inertial1(i,2),R_inertial1(i,3), attx(1), attx(2), attx(3), 'b')
    hold on
    quiver3(R_inertial1(i,1),R_inertial1(i,2),R_inertial1(i,3), atty(1), atty(2), atty(3), 'b')
    hold on
    quiver3(R_inertial1(i,1),R_inertial1(i,2),R_inertial1(i,3), attz(1), attz(2), attz(3), 'r')    
    xlabel('X axis'), ylabel('Y axis'), zlabel('Z axis')
    hold off
    legend('Trajectory', 'x_b', 'y_b', 'z_b','interpreter', 'TeX')
    grid on
    % pause(.00000001)
    xlim([min(R_inertial(:,1))-3 max(R_inertial(:,1))+3])
    ylim([min(R_inertial(:,2))-3 max(R_inertial(:,2))+3])
    zlim([min(R_inertial(:,3))-3 max(R_inertial(:,3))+3])
end

%% Attitude of the non-spinning disk w/ Tn
phi = sim.EulAng(:,1);
theta = sim.EulAng(:,2);
psi = sim.EulAng(:,3);
PHI = zeros(size(phi));
THETA = zeros(size(phi));
PSI = zeros(size(phi));
lambda = sim.lambda;
R_inertial = sim.R_inertial(:,:)';
% % Define the filename for the video
% videoFilename = 'trajectory_video_nospindisk.avi';

% % Create a VideoWriter object
% video = VideoWriter(videoFilename);
% 
% % Set the frame rate (frames per second)
% numFrames = length(phi(1:50:end));
% t_video = 10;                 % Desired video duration [s]
% frameRate = numFrames/t_video; % Adjust as needed
% video.FrameRate = frameRate;
% 
% % Open the VideoWriter
% open(video);

x_v = [1; 0; 0];
y_v = [0; 1; 0];
z_v = [0; 0; 1];
% R_inertial = zeros(size(R_inertial));
for i = 1:length(sim.tout)
    R1 = [1, 0, 0; 0, cos(phi(i)), sin(phi(i)); 0, -sin(phi(i)), cos(phi(i))];
    R2 = [cos(theta(i)), 0, -sin(theta(i)); 0, 1, 0; sin(theta(i)), 0, cos(theta(i))];
    R3 = [cos(psi(i)), sin(psi(i)), 0; -sin(psi(i)), cos(psi(i)), 0; 0, 0, 1];
    T0 = R1*R2*R3;
    Tn = [cos(lambda(i)) -sin(lambda(i)) 0; sin(lambda(i)) cos(lambda(i)) 0; 0 0 1];
    TI = Tn*T0*TI0;
    THETA(i) = -asin(TI(1,3));
    PSI(i) = acos( TI(1,1)/cos(THETA(i)) );
    PHI(i) = acos( TI(3,3)/cos(THETA(i)) );
end

figure(2)
plot(sim.tout, rad2deg(PHI),'LineWidth', 1)
xlabel('Time [s]')
ylabel('\Phi [deg]')
hold on
plot(sim.tout, rad2deg(THETA),'LineWidth', 1)
xlabel('Time [s]')
ylabel('Angle [deg]')
legend('\Phi', '\Theta')
figure(4)
plot(sim.tout, rad2deg(PSI),'LineWidth', 1)
xlabel('Time [s]')
ylabel('\Psi [deg]')
max(rad2deg(PHI))
%% Animation
for i = 1:50:length(sim.tout)
    R1 = [1, 0, 0; 0, cos(phi(i)), sin(phi(i)); 0, -sin(phi(i)), cos(phi(i))];
    R2 = [cos(theta(i)), 0, -sin(theta(i)); 0, 1, 0; sin(theta(i)), 0, cos(theta(i))];
    R3 = [cos(psi(i)), sin(psi(i)), 0; -sin(psi(i)), cos(psi(i)), 0; 0, 0, 1];
    T0 = R1*R2*R3;
    invT0 = transpose(T0);
    attx = invTI0*invT0*x_v;
    atty = invTI0*invT0*y_v;
    attz = invTI0*invT0*z_v;
    
    Tn = [cos(lambda(i)) -sin(lambda(i)) 0; sin(lambda(i)) cos(lambda(i)) 0; 0 0 1];
    invTn = transpose(Tn);
    attxn = invTI0*invT0*invTn*x_v;
    attyn = invTI0*invT0*invTn*y_v;
    attzn = invTI0*invT0*invTn*z_v;
    TI = Tn*T0*TI0;
    THETA(i) = -asin(TI(1,3));
    PSI(i) = acos( TI(1,1)/cos(THETA(i)) );
    PHI(i) = acos( TI(3,3)/cos(THETA(i)) );

    figure(1);
    plot3(R_inertial(:,1),R_inertial(:,2),R_inertial(:,3),'k-')
    hold on
    quiver3(R_inertial(i,1),R_inertial(i,2),R_inertial(i,3), attx(1), attx(2), attx(3), 'g')
    hold on
    quiver3(R_inertial(i,1),R_inertial(i,2),R_inertial(i,3), atty(1), atty(2), atty(3), 'g')
    hold on
    quiver3(R_inertial(i,1),R_inertial(i,2),R_inertial(i,3), attz(1), attz(2), attz(3), 'r')
    hold on
    quiver3(R_inertial(i,1),R_inertial(i,2),R_inertial(i,3), attxn(1), attxn(2), attxn(3), 'b')
    hold on
    quiver3(R_inertial(i,1),R_inertial(i,2),R_inertial(i,3), attyn(1), attyn(2), attyn(3), 'b')
    hold on
    quiver3(R_inertial(i,1),R_inertial(i,2),R_inertial(i,3), attzn(1), attzn(2), attzn(3), 'm')    
    xlabel('X axis'), ylabel('Y axis'), zlabel('Z axis')
    hold off
    legend('Trajectory', 'x_b', 'y_b', 'z_b', 'x_{nospin}', 'y_{nospin}', '','interpreter', 'TeX')
    grid on

    xlim([min(R_inertial(:,1))-1 max(R_inertial(:,1))+1])
    ylim([min(R_inertial(:,2))-1 max(R_inertial(:,2))+1])
    zlim([min(R_inertial(:,3))-1 max(R_inertial(:,3))+1])
    %  % Capture the current frame
    % frame = getframe(gcf);
    % 
    % % Write the frame to the video
    % writeVideo(video, frame);
end
% % Close the VideoWriter
% close(video);

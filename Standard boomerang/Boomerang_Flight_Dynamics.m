%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                   Boomerang trajectory simulation                   %%%
%%%                   Davide Di Santis - November 2023                  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Source: Flight Dynamics of the Boomerang, Part 1: Fundamental Analysis%%
%% Data
clear
close all
clc

tic
%--------------------------Simulation parameters--------------------------%
t = 5;                     % Simulation time [s]
dt = 2e-4;                 % Timestep for fixed-step simulation [s]

%---------------------------Planetary parameters--------------------------%
g = 9.81;               % Gravity acceleration [m/s^2]
rho = 1.22;             % Air density at sea level [kg/m^3]
v_wind = [0; 0; 0];  % Wind velocity [m/s]

%-----------------------Boomerang design parameters-----------------------%
l_blade = 3e-1; % Length of one blade [m]
nw = 2;         % Number of wings [ ]
R = 2.69e-1;    % Radius of rotation [m]
S = 2.28e-1;    % Disk area [m^2]
m_b = 1.30e-1;    % Boomerang mass [kg]
c = 4.88e-2;    % Mean chord [m]
m_add = zeros(3);
% m_add(3,3) = 8/3*rho*R^3;
m = m_b*eye(3) + m_add;
invm = inv(m);

x_ac = 6.1e-2; % Position of aerodynamic center in body coordinates
LAMBDAj = 120;  % Wing sweep angle [deg]
LAMBDAj = deg2rad(LAMBDAj);
gamma = 120;    % Folding angle [deg]
gamma = deg2rad(gamma);
betaj = 0;      % Wing coning angle 
betaj = deg2rad(betaj);
thetaj = 0;     % Wing pitch angle
thetaj = deg2rad(thetaj);

%------------------Moments of inertia of a single blade-------------------%
I_xi = 1.88e-3; 
I_eta = 4.88e-6;
I_zeta = 1.92e-3;
I_xieta = 5.14e-21;

Jj = [
    I_xi I_xieta 0;
    I_xieta I_eta 0;
    0 0 I_zeta
    ];

%----------------------------Launch conditions----------------------------%
U0 = [25*cosd(30); 25*sind(30); 0];  % Initial throw speed in inertial frame [m/s]
omega0 = [0; 0; 10*pi*2];            % Initial angular speed [rad/s]
PHI0 = 73;                           % Initial roll angle (Between non-spinning frame and inertial) [deg]
PHI0 = deg2rad(PHI0);
THETA0 = 0;                          % Initial pitch angle (Between non-spinning frame and inertial) [deg]
THETA0 = deg2rad(THETA0);
PSI0 = 225;                          % Initial yaw angle (Between non-spinning frame and inertial) [deg]
PSI0 = deg2rad(PSI0);
R_pos0 = [0; 0; 1.5];                % Initial position in inertial frame

%-------------Parameters for integration of aerodynamic forces------------%
n = 50; % Number of intervals
n = n+1;% Number of points
l_integrate = linspace(0,l_blade,n);

%---------------------------Rotation matrices-----------------------------%
Tj = zeros(3,3, nw);
invTj = zeros(3,3,nw);

for i = 0:(nw-1)
Rj1 = [1, 0, 0; 0, cos(betaj), sin(betaj); 0, -sin(betaj), cos(betaj)];
    Rj2 = [cos(thetaj), 0, -sin(thetaj);0, 1, 0; sin(thetaj), 0, cos(thetaj)];
    Rj3 = [cos(LAMBDAj-pi/2 + gamma*i), sin(LAMBDAj-pi/2 + gamma*i), 0; ...
        -sin(LAMBDAj-pi/2 + gamma*i), cos(LAMBDAj-pi/2+ gamma*i), 0;0, 0, 1];
    Tj(:,:,i+1) = Rj2*Rj1*Rj3;
    invTj(:,:,i+1) = transpose(Tj(:,:,i+1));
end

%-------------Computation of moment of inertia wrt body frame-------------%
Ji = zeros(3,3,nw);
J = zeros(3);

for i = 1:nw
    Ji(:,:,i) = invTj(:,:,i)*Jj*Tj(:,:,i);
    J = J + (Ji(:,:,i));
end

J(2,2) = J(2,2) + m_b*x_ac^2;
J(3,3) = J(3,3) + m_b*x_ac^2;

J = diag(diag(J));

invJ = inv(J);

RI01 = [1, 0, 0; 0, cos(PHI0), sin(PHI0); 0, -sin(PHI0), cos(PHI0)];
RI02 = [cos(THETA0), 0, -sin(THETA0); 0, 1, 0; sin(THETA0), 0, cos(THETA0)];
RI03 = [cos(PSI0), sin(PSI0), 0; -sin(PSI0), cos(PSI0), 0; 0, 0, 1];
TI0 = RI01*RI02*RI03;
invTI0 = transpose(TI0);

%--------------------Starting conditions, body frame----------------------%
u0 = TI0*U0;
r0 = TI0*R_pos0;

Eul0 = [0; 0; 0];
Quat0 = [1; 0; 0; 0];
%------------------------Aerodynamic coefficients-------------------------%
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

% sim = sim("Boomerang_Simulink_quat_beta.slx");
sim = sim("Boomerang_Simulink.slx");

toc
%% Plot results 

% phi = sim.EulAng(:,1);
% theta = sim.EulAng(:,2);
% psi = sim.EulAng(:,3);
R_inertial = sim.R_inertial(:,:)';
U = sim.U(:,:)';
normU = vecnorm(U');
% R_inertial = R_inertial;
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

plot(sim.tout,sim.omega(:,3))
xlabel('Time [s]')
ylabel('omega [rad/s]')
title('Spin rate over time')
%% Trajectory and attitude - Euler Angles
phi = sim.EulAng(:,1);
theta = sim.EulAng(:,2);
psi = sim.EulAng(:,3);
R_inertial = sim.R_inertial(:,:)';
% % Define the filename for the video
% videoFilename = 'trajectory_video.avi';

% % Create a VideoWriter object
% video = VideoWriter(videoFilename);
% 
% % Set the frame rate (frames per second)
% numFrames = length(phi(1:50:end));
% frameRate = numFrames/(sim.tout(end)); % Adjust as needed
% video.FrameRate = frameRate;
% 
% % Open the VideoWriter
% open(video);

x_v = [1; 0; 0];
y_v = [0; 1; 0];
z_v = [0; 0; 1];
% R_inertial = zeros(size(R_inertial));
for i = 1:50:length(sim.tout)
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

    figure(1);
    plot3(R_inertial(:,1),R_inertial(:,2),R_inertial(:,3),'k-')
    hold on
    quiver3(R_inertial(i,1),R_inertial(i,2),R_inertial(i,3), attx(1), attx(2), attx(3), 'b')
    hold on
    quiver3(R_inertial(i,1),R_inertial(i,2),R_inertial(i,3), atty(1), atty(2), atty(3), 'b')
    hold on
    quiver3(R_inertial(i,1),R_inertial(i,2),R_inertial(i,3), attz(1), attz(2), attz(3), 'r')    
    xlabel('X axis'), ylabel('Y axis'), zlabel('Z axis')
    hold off
    legend('Trajectory', 'x_b', 'y_b', 'z_b','interpreter', 'TeX')
    title('Time elapsed: ', num2str(sim.tout(i)))
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

%% Trajectory and attitude - Quaternions
% phi = sim.EulAng(1,:);
% theta = sim.EulAng(2,:);
% psi = sim.EulAng(3,:);
R_inertial = sim.R_inertial(:,:)';
% % Define the filename for the video
% videoFilename = 'trajectory_video.avi';

% % Create a VideoWriter object
% video = VideoWriter(videoFilename);
% 
% % Set the frame rate (frames per second)
% numFrames = length(phi(1:50:end));
% frameRate = numFrames/(sim.tout(end)); % Adjust as needed
% video.FrameRate = frameRate;
% 
% % Open the VideoWriter
% open(video);

x_v = [1; 0; 0];
y_v = [0; 1; 0];
z_v = [0; 0; 1];
for i = 1:100:length(sim.tout)
    % R1 = [
    % 1, 0, 0;
    % 0, cos(phi(i)), sin(phi(i));
    % 0, -sin(phi(i)), cos(phi(i))
    % ];
    % R2 = [
    % cos(theta(i)), 0, -sin(theta(i));
    % 0, 1, 0;
    % sin(theta(i)), 0, cos(theta(i))
    % ];
    % R3 = [
    % cos(psi(i)), sin(psi(i)), 0;
    % -sin(psi(i)), cos(psi(i)), 0;
    % 0, 0, 1
    % ];
    % T0 = R1*R2*R3;
    % invT0 = transpose(T0);
    invT0 = transpose(sim.DCM(:,:,i));
    attx = invTI0*invT0*x_v;
    atty = invTI0*invT0*y_v;
    attz = invTI0*invT0*z_v;
    
    figure(1);
    plot3(R_inertial(:,1),R_inertial(:,2),R_inertial(:,3),'k-')
    hold on
    quiver3(R_inertial(i,1),R_inertial(i,2),R_inertial(i,3), attx(1), attx(2), attx(3), 'b')
    hold on
    quiver3(R_inertial(i,1),R_inertial(i,2),R_inertial(i,3), atty(1), atty(2), atty(3), 'b')
    hold on
    quiver3(R_inertial(i,1),R_inertial(i,2),R_inertial(i,3), attz(1), attz(2), attz(3), 'r')    
    xlabel('X axis'), ylabel('Y axis'), zlabel('Z axis')
    hold off
    legend('Trajectory', 'x_b', 'y_b', 'z_b','interpreter', 'TeX')
    grid on
    % pause(.00000001)
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

%% Trajectory and attitude of the non-spinning disk w/ Tn
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

%% Compare IMU Data
close all

R_inertial = sim.R_inertial(:,:)';
U = sim.U(:,:)';
normU = vecnorm(U');
u = sim.u(:,:);
acc = sim.acc(:,:);
sample = 1/(f*dt);
t_s = linspace(0,t(end), ceil(f*sim.tout(end)));

R_s = R_inertial(1:sample:end,:);
u_s = u(:,1:sample:end);
acc_s = acc(:,1:sample:end);

acc_nav = sim.acc_nav;
u_nav = sim.u_nav;
R_nav = sim.R_nav(:,:);

err_vel = vecnorm(u_s - u_nav');
err_pos = vecnorm(R_s' - R_nav);

figure(1)
plot(sim.tout, acc(1,:), 'k-')
hold on
plot(t_s, acc_nav(:,1), 'k--')
legend('Dynamic model data', 'IMU estimates')
xlabel('Time [s]')
ylabel('Acceleration [m/s^2]')
title('Acceleration along x - Model vs estimate comparison')

figure(2)
plot(sim.tout, acc(2,:), 'k-')
hold on
plot(t_s, acc_nav(:,2), 'k--')
legend('Dynamic model data', 'IMU estimates')
xlabel('Time [s]')
ylabel('Acceleration [m/s^2]')
title('Acceleration along y - Model vs estimate comparison')

figure(3)
plot(sim.tout, acc(3,:), 'k-')
hold on
plot(t_s, acc_nav(:,3), 'k--')
legend('Dynamic model data', 'IMU estimates')
xlabel('Time [s]')
ylabel('Acceleration [m/s^2]')
title('Acceleration along z - Model vs estimate comparison')

figure(4)
plot(sim.tout, normU, 'k-')
hold on
plot(t_s, vecnorm(u_nav'), 'k--')
legend('Dynamic model data', 'IMU estimates')
xlabel('Time [s]')
ylabel('Speed [m/s]')
title('Speed - Model vs estimate comparison')

figure(5)
plot(R_inertial(:,1),R_inertial(:,2), 'k-')
hold on
plot(R_nav(1,:), R_nav(2,:), 'k--')
legend('Dynamic model data', 'IMU estimates')
axis equal
xlabel('Displacemenent along x')
ylabel('Displacemenent along y')
title('Top view trajectory - Model vs estimate comparison')

% figure(6)
% plot(sim.tout, err_vel)

figure(7)
plot(t_s, err_pos, 'k-')
xlabel('Time [s]')
ylabel('Error [m]')
title('Error in position estimation')

%% Compare spinning and non-spinning disk attitudes
phi = sim.EulAng(:,1);
theta = sim.EulAng(:,2);
psi = sim.EulAng(:,3);
PHI = sim.EULANG(1,:);
THETA = sim.EULANG(2,:);
PSI = sim.EULANG(3,:);
R_inertial = sim.R_inertial(:,:)';
% % Define the filename for the video
% videoFilename = 'trajectory_video.avi';

% % Create a VideoWriter object
% video = VideoWriter(videoFilename);
% 
% % Set the frame rate (frames per second)
% numFrames = length(phi(1:50:end));
% frameRate = numFrames/(sim.tout(end)); % Adjust as needed
% video.FrameRate = frameRate;
% 
% % Open the VideoWriter
% open(video);

x_v = [1; 0; 0];
y_v = [0; 1; 0];
z_v = [0; 0; 1];

for i = 1:50:length(sim.tout)
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
    
    R1n = [
    1, 0, 0;
    0, cos(PHI(i)), sin(PHI(i));
    0, -sin(PHI(i)), cos(PHI(i))
    ];
    R2n = [
    cos(THETA(i)), 0, -sin(THETA(i));
    0, 1, 0;
    sin(THETA(i)), 0, cos(THETA(i))
    ];
    R3n = [
    cos(PSI(i)), sin(PSI(i)), 0;
    -sin(PSI(i)), cos(PSI(i)), 0;
    0, 0, 1
    ];
    T = R1n*R2n*R3n;
    invT = transpose(T);
    
    attx = invTI0*invT0*x_v;
    atty = invTI0*invT0*y_v;
    attz = invTI0*invT0*z_v;
    attxn = invT*x_v;
    attyn = invT*y_v;
    attzn = invT*z_v;

    figure(1);
    plot3(R_inertial(:,1),R_inertial(:,2),R_inertial(:,3),'k-')
    hold on
    quiver3(R_inertial(i,1),R_inertial(i,2),R_inertial(i,3), attx(1), attx(2), attx(3), 'b')
    hold on
    quiver3(R_inertial(i,1),R_inertial(i,2),R_inertial(i,3), atty(1), atty(2), atty(3), 'b')
    hold on
    quiver3(R_inertial(i,1),R_inertial(i,2),R_inertial(i,3), attz(1), attz(2), attz(3), 'r')    
    hold on
    quiver3(R_inertial(i,1),R_inertial(i,2),R_inertial(i,3), attxn(1), attxn(2), attxn(3), 'c')
    hold on
    quiver3(R_inertial(i,1),R_inertial(i,2),R_inertial(i,3), attyn(1), attyn(2), attyn(3), 'c')
    hold on
    quiver3(R_inertial(i,1),R_inertial(i,2),R_inertial(i,3), attzn(1), attzn(2), attzn(3), 'g')    
    xlabel('X axis'), ylabel('Y axis'), zlabel('Z axis')
    hold off
    legend('Trajectory', 'x_b', 'y_b', 'z_b','interpreter', 'TeX')
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
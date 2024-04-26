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
m = 1.30e-1;    % Boomerang mass [kg]
c = 4.88e-2;    % Mean chord [m]

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
PSI0 = 225;                           % Initial yaw angle (Between non-spinning frame and inertial) [deg]
PSI0 = deg2rad(PSI0);
THETA0 = 0;                          % Initial pitch angle (Between non-spinning frame and inertial) [deg]
THETA0 = deg2rad(THETA0);  
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

J(2,2) = J(2,2) + m*x_ac^2;
J(3,3) = J(3,3) + m*x_ac^2;

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

%---------------------Desired trajectory: returning path------------------% 
load('Returning_Traj_XYZ-U-tout.mat');
x_des = R_inertial1(:,1);
y_des = R_inertial1(:,2);

X_des = [x_des, y_des];

% %-------------------Desired trajectory: circular------------------------%
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

%---------------------Mechanical control parameters-----------------------%
servo_rpm = 500;           % Rotation speed of the servo/rotation speed of the wing
servo_rpm = deg2rad(servo_rpm);
theta_sat = 25;            % Max/min rotation of the wing [deg]
theta_sat = deg2rad(theta_sat);

%-----------------------------PID parameters------------------------------%
Kp = 0; 
Ki = 0;
Kd = 0;
Kp_roll = 2;
Ki_roll = 15;
Kd_roll = .004;
PHI_des = 73;    % Desired angle between the horizon plane and z_body axis
PHI_des = deg2rad(PHI_des);
sim = sim("Boomerang_pitch_att_simulink.slx");

toc
%% Plot top view trajectory only
% u = zeros(length(sim.tout),3);
% normU = zeros(length(sim.tout),1);

%-Load the correct file for uncontrolled trajectory under effects of wind-%
load('Wind_wx_Traj_XYZ-U-tout.mat')

x_w = R_inertial_w(:,1);
y_w = R_inertial_w(:,2);

R_inertial = sim.R_inertial(:,:);
% u(i,:) = sim.U(:,:);
% normU = vecnorm(u);

figure(1)
plot(R_inertial(1,1),R_inertial(2,1), 'go', ...
    R_inertial(1,:), R_inertial(2,:), 'k-', ...
    R_inertial(1,end), R_inertial(2,end), 'ro')
xlabel('Displacement along x [m]')
ylabel('Displacement along y [m]')
title('Trajectory - top view (xy plane)')
hold on
plot(x_des,y_des, 'b-');
hold on
plot(x_w, y_w, 'r--');
legend('', 'PD trajectory', '', 'Desired trajectory','Uncontrolled trajectory', 'Location','southeast')
axis equal
%% Trajectory and attitude
phi = sim.EulAng(:,1);
theta = sim.EulAng(:,2);
psi = sim.EulAng(:,3);
R_inertial = sim.R_inertial(:,:)';

x_v = [1; 0; 0];
y_v = [0; 1; 0];
z_v = [0; 0; 1];
% R_inertial = zeros(size(R_inertial));
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
    
    figure(1);
    plot3(R_inertial(:,1),R_inertial(:,2),R_inertial(:,3),'k-')
    hold on
    quiver3(R_inertial(i,1),R_inertial(i,2),R_inertial(i,3), attx(1), attx(2), attx(3), 'b')
    hold on
    quiver3(R_inertial(i,1),R_inertial(i,2),R_inertial(i,3), atty(1), atty(2), atty(3), 'b')
    hold on
    quiver3(R_inertial(i,1),R_inertial(i,2),R_inertial(i,3), attz(1), attz(2), attz(3), 'r')    
    xlabel('X axis'), ylabel('Y axis'), zlabel('Z axis')
    hold on
    quiver3(R_inertial(i,1),R_inertial(i,2),R_inertial(i,3), 0, 0, 1, 'g') 
    hold off
    legend('Trajectory', 'x_b', 'y_b', 'z_b', 'Z_{des}','interpreter', 'TeX', 'FontSize', 11)
    grid on
    % pause(.00000001)
    xlim([min(R_inertial(:,1))-1 max(R_inertial(:,1))+1])
    ylim([min(R_inertial(:,2))-1 max(R_inertial(:,2))+1])
    zlim([min(R_inertial(:,3))-1 max(R_inertial(:,3))+1])
end
%% Plot roll angle and desired roll angle
close all
clc
phi = sim.EulAng(:,1);
theta = sim.EulAng(:,2);
psi = sim.EulAng(:,3);
THETA = zeros(size(phi));
PHI = zeros(size(phi));
lambda = sim.lambda;
load('PHI_natural')
for i = 1:length(sim.tout)
    R1 = [1, 0, 0; 0, cos(phi(i)), sin(phi(i)); 0, -sin(phi(i)), cos(phi(i))];
    R2 = [cos(theta(i)), 0, -sin(theta(i)); 0, 1, 0; sin(theta(i)), 0, cos(theta(i))];
    R3 = [cos(psi(i)), sin(psi(i)), 0; -sin(psi(i)), cos(psi(i)), 0; 0, 0, 1];
    T0 = R1*R2*R3;
    Tn = [cos(lambda(i)) -sin(lambda(i)) 0; sin(lambda(i)) cos(lambda(i)) 0; 0 0 1];
    TI = Tn*T0*TI0;
    THETA(i) = -asin(TI(1,3));
    % PSI(i) = acos( TI(1,1)/cos(THETA(i)) );
    PHI(i) = acos( TI(3,3)/cos(THETA(i)) );
end
% figure(2)
plot(sim.tout, rad2deg(PHI), sim.tout, rad2deg(PHI_des)*ones(size(sim.tout)),'LineWidth', 1)
hold on
plot(t1, rad2deg(PHI1), 'k--')
legend('Roll angle', 'Desired roll angle', 'Uncontrolled roll angle')
% figure(3)
% plot(sim.tout, rad2deg(THETA),'LineWidth', 1)
% figure(4)
% plot(sim.tout, rad2deg(PSI),'LineWidth', 1)

% Averages
[~,maxind] = findpeaks(lambda); % Find the local maxima of lambda 
maxind = [1;maxind];
mean_t = zeros(length(maxind)-1,1);
mean_PHI = zeros(length(maxind)-1,1);

for i = 1:(length(maxind)-1)
    mean_t(i) = mean(sim.tout(maxind(i):maxind(i+1)));
    mean_PHI(i) = mean(PHI(maxind(i):maxind(i+1)));
end

figure(2)
plot(sim.tout,rad2deg(PHI_des)*ones(size(sim.tout)), mean_t, rad2deg(mean_PHI), 'LineWidth',1.5);
hold on
plot(mean_t1, rad2deg(mean_PHI1), 'k--','LineWidth', 0.5)
xlim([0, sim.tout(end)])
legend('Desired angle', 'Avg Actual angle', 'Avg Uncontrolled angle - natural trajectory', 'FontSize', 11)
xlabel('Time [s]')
ylabel('Roll angle \Phi [deg]')

s = abs(max(abs(mean_PHI)) - PHI_des);
fprintf('\nThe maximum overshoot is %f degrees \n', rad2deg(s))

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
%% Plot all trajectory results 
R_inertial = sim.R_inertial(:,:);
u = sim.U(:,:);
normU = vecnorm(u);

% Plots
figure(1)
plot(sim.tout, R_inertial(3,:), 'k-')
xlabel('Time [s]')
ylabel('Altitude [m]')

figure(2)
plot(R_inertial(1,1),R_inertial(2,1), 'go', ...
    R_inertial(1,:), R_inertial(2,:), 'k-', ...
    R_inertial(1,end), R_inertial(2,end), 'ro')
xlabel('Displacement along x [m]')
ylabel('Displacement along y [m]')
title('Trajectory - top view (xy plane)')
axis equal
% xlim([-10 10])
% ylim([-5 18])

figure(3)
plot(R_inertial(1,1),R_inertial(3,1), 'go', ...
    R_inertial(1,:), R_inertial(3,:), 'k-', ...
    R_inertial(1,end), R_inertial(3,end), 'ro')
xlabel('Displacement along x [m]')
ylabel('Displacement along z [m]')
title('Trajectory - xz plane')
axis equal
% xlim([-10 10])
% ylim([-2 8])

figure(4)
plot(R_inertial(2,1), R_inertial(3,1), 'go', ...
    R_inertial(2,:), R_inertial(3,:), 'k-', ...
    R_inertial(2,end), R_inertial(3,end), 'ro')
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

figure(6)
plot3(R_inertial(1,1), R_inertial(2,1), R_inertial(3,1), 'go', ...
    R_inertial(1,:), R_inertial(2,:), R_inertial(3,:), 'k-', ...
    R_inertial(1,end), R_inertial(2,end), R_inertial(3,end), 'ro')
title('3D Trajectory')
grid on

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
% videoFilename = 'trajectory_video_lambda.avi';

% Create a VideoWriter object
% video = VideoWriter(videoFilename);

% % Set the frame rate (frames per second)
% numFrames = length(phi(1:50:end));
% t_video = 20;                 % Desired video duration [s]
% frameRate = numFrames/t_video; % Adjust as needed
% video.FrameRate = frameRate;

% Open the VideoWriter
% open(video);

x_v = [1; 0; 0];
y_v = [0; 1; 0];
z_v = [0; 0; 1];
% R_inertial = zeros(size(R_inertial));


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

    figure(1);
    plot3(R_inertial(:,1),R_inertial(:,2),R_inertial(:,3),'k-')
    hold on
    quiver3(R_inertial(i,1),R_inertial(i,2),R_inertial(i,3), attx(1), attx(2), attx(3), 'b')
    hold on
    quiver3(R_inertial(i,1),R_inertial(i,2),R_inertial(i,3), atty(1), atty(2), atty(3), 'b')
    hold on
    quiver3(R_inertial(i,1),R_inertial(i,2),R_inertial(i,3), attz(1), attz(2), attz(3), 'r')
    hold on
    quiver3(R_inertial(i,1),R_inertial(i,2),R_inertial(i,3), attxn(1), attxn(2), attxn(3), 'g')
    hold on
    quiver3(R_inertial(i,1),R_inertial(i,2),R_inertial(i,3), attyn(1), attyn(2), attyn(3), 'g')
    hold on
    quiver3(R_inertial(i,1),R_inertial(i,2),R_inertial(i,3), attzn(1), attzn(2), attzn(3), 'm')    
    xlabel('X axis'), ylabel('Y axis'), zlabel('Z axis')
    hold off
    legend('Trajectory', 'x_b', 'y_b', 'z_b','interpreter', 'TeX')
    grid on

    xlim([min(R_inertial(:,1))-1 max(R_inertial(:,1))+1])
    ylim([min(R_inertial(:,2))-1 max(R_inertial(:,2))+1])
    zlim([min(R_inertial(:,3))-1 max(R_inertial(:,3))+1])
     % Capture the current frame
    % frame = getframe(gcf);
    % 
    % % Write the frame to the video
    % writeVideo(video, frame);
end
% % Close the VideoWriter
% close(video);
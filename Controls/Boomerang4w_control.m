%% Data
clear
% close all
clc

tic
%--------------------------Simulation parameters--------------------------%
t = 2.85;                     % Simulation time [s]
dt = 2e-4;                    % Fixed-step timestep
f_con = 1000;                 % Control frequency [Hz]                
dt_con = 1/f_con;             % Control timestep [s]
minstep = 1e-5;               % Min timestep for variable-step simulation [s]
maxstep = 1.66e-4;            % Max timestep for variable-step simulation [s]

%---------------------------Planetary parameters--------------------------%
g = 9.81;                % Gravity acceleration [m/s^2]
rho = 1.22;              % Air density at sea level [kg/m^3]
v_wind = [-1; 0; 0];      % Wind speed in inertial reference frame [m/s]

%-----------------------Boomerang design parameters-----------------------%
l_blade = 25e-2;       % Length of one blade [m]
nw = 4;                 % Number of wings [ ]
R = l_blade;            % Radius of rotation [m]
S = pi*R^2;             % Disk area [m^2]
m = 200e-3;          % Boomerang mass [kg]
c = 5e-2;               % Mean chord [m]

d_ac = c/4; % Position of aerodynamic center in body coordinates
LAMBDAj = 0;  % Wing sweep angle [deg]
LAMBDAj = deg2rad(LAMBDAj);
gamma = 90;    % Folding angle [deg]
gamma = deg2rad(gamma);
betaj = 4;      % Wing coning angle 
betaj = deg2rad(betaj);
thetaj = 25;     % Wing pitch angle at root
x_ac = zeros(3,nw);
for i = 1:nw
    Rac = [cos(gamma*i), sin(gamma*i), 0; -sin(gamma*i), cos(gamma*i), 0;0, 0, 1];
    x_ac(:,i) = Rac*[0; -d_ac; 0];
end

%-----------------Moments of inertia of a single blade--------------------%
I_xi = (m/nw)/12*(l_blade^2 + (0.1*c)^2) + m/nw*(l_blade/2)^2; 
I_eta = (m/nw)/12*(c^2);
I_zeta = I_xi+I_eta;%(m_b/nw)/12*(l_blade^2 + c^2) + m_b/nw*(l_blade/2)^2;
I_xieta = 0;

Jj = [
    I_xi I_xieta 0;
    I_xieta I_eta 0;
    0 0 I_zeta
    ];

%-------------------------Launch conditions-------------------------------%
ThAng = 45;                                 % Throw angle from East direction
U0mod = 25;
U0 = [U0mod*cosd(ThAng); U0mod*sind(ThAng); 0];   % Initial throw speed in inertial frame [m/s]
omega0 = [0; 0; 10*pi*2];                   % Initial angular speed [rad/s]
PHI0 = 75;                                  % Initial roll angle (Between non-spinning frame and inertial) [deg]
PHI0 = deg2rad(PHI0);
THETA0 = 0;                                 % Initial pitch angle (Between non-spinning frame and inertial) [deg]
THETA0 = deg2rad(THETA0);  
PSI0 = 220;                                 % Initial yaw angle (Between non-spinning frame and inertial) [deg]
PSI0 = deg2rad(PSI0);
R_pos0 = [0; 0; 1.5];                       % Initial position in inertial frame

%------------Parameters for integration of aerodynamic forces-------------%
n = 50; % Number of intervals
n = n + 1; 
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

%------------Computation of moment of inertia wrt body frame--------------%
Ji = zeros(3,3,nw);
J = zeros(3);

for i = 1:nw
    Ji(:,:,i) = invTj(:,:,i)*Jj*Tj(:,:,i);
    J = J + (Ji(:,:,i));
    % for j = 1:3
    %     for k = 1:3
    %         J(j,k) = J(j,k) + m_b/nw*(norm(x_ac(:,nw))^2 - x_ac(j,nw)*x_ac(k,nw));
    %     end
    % end
end

% J(2,2) = J(2,2) + m_b*x_ac^2;
% J(3,3) = J(3,3) + m_b*x_ac^2;

% J = diag(diag(J));

invJ = inv(J);

%-------------------Starting conditions in body frame---------------------%
Eul0 = [0; 0; 0];
Quat0 = [1; 0; 0; 0];
u0 = TI0*U0;
r0 = TI0*R_pos0;

%------------------------Aerodynamic coefficients-------------------------%
% CL data
M1 = readmatrix("Cl.csv");
x_CL = M1(:,3);
y_CL = 1.35*M1(:,2);
% y_CL = M1(:,2);
CLdata = [x_CL, y_CL];

% CD data
M2 = readmatrix("Cd.csv");
x_CD = M2(:,3);
y_CD = 1.35*M2(:,2);
% y_CD = M2(:,2);
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

%---------------------Desired trajectory: returning path------------------% 
load('Returning_Traj_4wings.mat');
x_des = R_inertial1(:,1);
y_des = R_inertial1(:,2);

X_des = [x_des, y_des];
%-------------------------Pitch control parameters------------------------%
servo_rpm = 500;           % Rotation speed of the servo/rotation speed of the wing
servo_rpm = deg2rad(servo_rpm);
theta_sat = 35;            % Max/min rotation of the wing [deg]
theta_sat = deg2rad(theta_sat);
%--------------------------CoM control parameters-------------------------%
CoM_sat = 4e-2;            % Max/min shift of CoM [m]
rate = 1;                  % Speed at which the CoM can shift [m/s]
%-----------------------Dihedral control parameters-----------------------%
max_dd = 45;              % Max delta of joint angle [deg]
max_dd = deg2rad(max_dd);  
min_dd = -45;              % Min delta of joint angle [deg]
min_dd = deg2rad(min_dd);
% %-------------------------Pitch PID parameters----------------------------%
Kp = .0; 
Ki = 0.;
Kd = 0.;
% %-----------------------------CoM PD parameters---------------------------%
% Kp = -.02;
% Ki = 0;
% Kd = -.1;
%-----------------------------Dihedral PD parameters----------------------%
% Kp = 0.5;
% Ki = 0;
% Kd = 0.2;
% Kp = 0; Ki = 0; Kd = 0;
%---------------Cyclic pitch - attitude control PID parameters------------%
Kp_roll = 0;
Ki_roll = 0;
Kd_roll = 0;
PHI_des = 85;    % Desired angle between the horizon plane and z_body axis
PHI_des = deg2rad(PHI_des);

sim = sim("Boomerang_pitch_att_simulink.slx");
% sim = sim("Boomerang_CoM_simulink");
% sim = sim("Boomerang_dihedral_simulink");

toc
%% Plot top view trajectory and compare with desired and uncontrolled trajectories
% u = zeros(length(sim.tout),3);
% normU = zeros(length(sim.tout),1);

%-Load the correct file for uncontrolled trajectory under effects of wind-%
load('Wind_wx_4wings.mat')

x_w = R_inertial_w(1,:);
y_w = R_inertial_w(2,:);

R_inertial = sim.R_inertial(:,:);
fprintf('\nThe final distance from the thrower is %f m\n',norm(R_inertial(:,end)))

figure(1)
plot(R_inertial(1,1),R_inertial(2,1), 'go', ...
    R_inertial(1,:), R_inertial(2,:), 'k-', ...
    R_inertial(1,end), R_inertial(2,end), 'ro')
xlabel('Displacement along x [m]')
ylabel('Displacement along y [m]')
title('Trajectory - top view (xy plane)')
hold on
plot(x_des,y_des, 'b--');
hold on
plot(x_w, y_w, 'r--');
legend('', 'Controlled trajectory', '', 'Desired trajectory','Uncontrolled trajectory', 'Location','southeast','FontSize',11)
axis equal
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
% xlim([0 5])
% ylim([0 25])

% figure(6)
% plot(sim.tout, sim.omega(:,3),'k')

% figure(6)
% plot3(R_inertial(1,1), R_inertial(1,2), R_inertial(1,3), 'go', ...
%     R_inertial(:,1), R_inertial(:,2), R_inertial(:,3), 'k-', ...
%     R_inertial(end,1), R_inertial(end,2), R_inertial(end,3), 'ro')
% title('3D Trajectory')
% grid on
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

%% Trajectory and attitude
phi = sim.EulAng(:,1);
theta = sim.EulAng(:,2);
psi = sim.EulAng(:,3);
R_inertial = sim.R_inertial(:,:)';
PHI = sim.PHI;
figure(2)
plot(sim.tout,zeros(size(sim.tout)), sim.tout, rad2deg(PHI));
legend('Desired', 'Actual')
x_v = [1; 0; 0];
y_v = [0; 1; 0];
z_v = [0; 0; 1];
R_inertial = zeros(size(R_inertial));
for i = 1:10:length(phi)
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
    legend('Trajectory', 'x_b', 'y_b', 'z_b', 'Z_des','interpreter', 'TeX', 'FontSize', 11)
    grid on
    % pause(.00000001)
    xlim([min(R_inertial(:,1))-1 max(R_inertial(:,1))+1])
    ylim([min(R_inertial(:,2))-1 max(R_inertial(:,2))+1])
    zlim([min(R_inertial(:,3))-1 max(R_inertial(:,3))+1])
end
%% Plot and compare energies with undisturbed returning trajectory

figure(1)
plot(sim.P_E, 'k-')
hold on
plot(P_Enc + 0*m*g, 'k--') % Adding 3.5*m*g to shift the curve as starting at 5.3 m 
title('Potential Energy over time')
xlabel('Time [s]')
ylabel('Energy [J]')
legend('E_P control', 'E_P no control', 'FontSize', 11)

figure(2)
plot(sim.K_E, 'k-')
hold on
plot(K_Enc, 'k--')
title('Kinetic Energy over time')
xlabel('Time [s]')
ylabel('Energy [J]')
legend('E_K control', 'E_K no control', 'FontSize', 11)

figure(3)
plot(sim.R_E, 'k-')
hold on
plot(R_Enc, 'k--')
title('Rotational Energy over time')
xlabel('Time [s]')
ylabel('Energy [J]')
legend('E_R control', 'E_R no control', 'FontSize', 11)

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
figure(3)
plot(sim.tout, rad2deg(THETA),'LineWidth', 1)
figure(4)
plot(sim.tout, rad2deg(PSI),'LineWidth', 1)

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
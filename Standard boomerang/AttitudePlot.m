%% Initial condition
clear
close all
clc

x_v = [1; 0; 0];
y_v = [0; 1; 0];
z_v = [0; 0; 1];
origin = [0; 0; 0];

PHI0 = 60;                                  % Initial roll angle (Between non-spinning frame and inertial) [deg]
PHI0 = deg2rad(PHI0);
THETA0 = 0;                                 % Initial pitch angle (Between non-spinning frame and inertial) [deg]
THETA0 = deg2rad(THETA0);
PSI0 = 270;                                 % Initial yaw angle (Between non-spinning frame and inertial) [deg]
PSI0 = deg2rad(PSI0);
ThAng = deg2rad(90);
speed = [cos(ThAng) sin(ThAng) 0];
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

attx = invTI0*x_v;
atty = invTI0*y_v;
attz = invTI0*z_v;

figure(1)
quiver3(origin(1), origin(2), origin(3), attx(1), attx(2), attx(3), 'b')
hold on
quiver3(origin(1), origin(2), origin(3), atty(1), atty(2), atty(3), 'b')
hold on
quiver3(origin(1), origin(2), origin(3), attz(1), attz(2), attz(3), 'r')
hold on
quiver3(origin(1), origin(2), origin(3), speed(1), speed(2), speed(3))
grid on;
xlabel('X axis'), ylabel('Y axis'), zlabel('Z axis')
legend('x_b', 'y_b', 'z_b','speed')
%% Dynamic plot w/ trajectory
clear
clc
close all

PHI0 = 70;                                  % Initial roll angle (Between non-spinning frame and inertial) [deg]
PHI0 = deg2rad(PHI0);
PSI0 = 240;                                 % Initial yaw angle (Between non-spinning frame and inertial) [deg]
PSI0 = deg2rad(PSI0);
THETA0 = 0;                                 % Initial pitch angle (Between non-spinning frame and inertial) [deg]
THETA0 = deg2rad(THETA0);  
load("Returning_Traj_XYZ-U-tout.mat");

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
phi = phi1;
theta = theta1;
psi = psi1;
% figure(1)

x_v = [1; 0; 0];
y_v = [0; 1; 0];
z_v = [0; 0; 1];
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
    invT0 = R3*R2*R1;
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
    xlim([-15 15])
    ylim([-15 10])
    zlim([-2 10])
end
%% Dynamic plot no trajectory
clear
clc
close all

PHI0 = 70;                                  % Initial roll angle (Between non-spinning frame and inertial) [deg]
PHI0 = deg2rad(PHI0);
PSI0 = 240;                                 % Initial yaw angle (Between non-spinning frame and inertial) [deg]
PSI0 = deg2rad(PSI0);
THETA0 = 0;                                 % Initial pitch angle (Between non-spinning frame and inertial) [deg]
THETA0 = deg2rad(THETA0);  
% load("Returning_Traj_XYZ-U-tout.mat");

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
% phi = phi1;
% theta = theta1;
% psi = psi1;
% figure(1)

x_v = [1; 0; 0];
y_v = [0; 1; 0];
z_v = [0; 0; 1];
origin = [0; 0; 0];
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
    invT0 = R3*R2*R1;
    attx = invTI0*invT0*x_v;
    atty = invTI0*invT0*y_v;
    attz = invTI0*invT0*z_v;
    figure(1)
    quiver3(origin(1),origin(2),origin(3), attx(1), attx(2), attx(3), 'b')
    hold on
    quiver3(origin(1),origin(2),origin(3), atty(1), atty(2), atty(3), 'b')
    hold on
    quiver3(origin(1),origin(2),origin(3), attz(1), attz(2), attz(3), 'r')
    
    hold off
    legend('x_b', 'y_b', 'z_b','interpreter', 'TeX')
    grid on
    % pause(.00000001)
    xlim([-2 2])
    ylim([-2 2])
    zlim([-2 2])
end

%% Dynamic plot no trajectory
clear
clc
close all

PHI0 = 70;                                  % Initial roll angle (Between non-spinning frame and inertial) [deg]
PHI0 = deg2rad(PHI0);
PSI0 = 240;                                 % Initial yaw angle (Between non-spinning frame and inertial) [deg]
PSI0 = deg2rad(PSI0);
THETA0 = 0;                                 % Initial pitch angle (Between non-spinning frame and inertial) [deg]
THETA0 = deg2rad(THETA0);  
load("Returning_Traj_XYZ-U-tout.mat");

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
phi = phi1;
theta = theta1;
psi = psi1;
% figure(1)

x_v = [1; 0; 0];
y_v = [0; 1; 0];
z_v = [0; 0; 1];
origin = [0; 0; 0];
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
    attx = (T0*TI0)*x_v;
    atty = (T0*TI0)*y_v;
    attz = (T0*TI0)*z_v;
    figure(1)
    quiver3(origin(1),origin(2),origin(3), attx(1), attx(2), attx(3), 'b')
    hold on
    quiver3(origin(1),origin(2),origin(3), atty(1), atty(2), atty(3), 'b')
    hold on
    quiver3(origin(1),origin(2),origin(3), attz(1), attz(2), attz(3), 'r')
    
    hold off
    legend('x_b', 'y_b', 'z_b','interpreter', 'TeX')
    grid on
    % pause(.00000001)
    xlim([-2 2])
    ylim([-2 2])
    zlim([-2 2])
end
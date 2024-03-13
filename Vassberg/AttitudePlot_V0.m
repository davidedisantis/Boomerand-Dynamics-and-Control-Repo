clear
clc
close all
x_v = [1; 0; 0];
y_v = [0; 1; 0];
z_v = [0; 0; 1];
origin = [0; 0; 0];

PHI0 = 70;                                  % Initial roll angle (Between non-spinning frame and inertial) [deg]
PHI0 = deg2rad(PHI0);
PSI0 = 240;                                 % Initial yaw angle (Between non-spinning frame and inertial) [deg]
PSI0 = deg2rad(PSI0);
THETA0 = 0;                                 % Initial pitch angle (Between non-spinning frame and inertial) [deg]
THETA0 = deg2rad(THETA0);  
load("Vassberg");

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
T = RI01*RI02*RI03;
invTI0 = transpose(T);

figure()

for i = 1:length(t(cut))
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
    attx = invT0*invTI0*x_v;
    atty = invT0*invTI0*y_v;
    attz = invT0*invTI0*z_v;
    figure(1)
    plot3([origin(1) attx(1)],[origin(2) attx(2)],[origin(3) attx(3)],'r-^', 'LineWidth',3);
    hold on
    plot3([origin(1) atty(1)],[origin(2) atty(2)],[origin(3) atty(3)],'r-^', 'LineWidth',3);
    hold on
    plot3([origin(1) attz(1)],[origin(2) attz(2)],[origin(3) attz(3)],'b-^', 'LineWidth',3);
    grid on;
    xlabel('X axis'), ylabel('Y axis'), zlabel('Z axis')
    legend('x_b', 'y_b', 'z_b')
    set(gca,'CameraPosition',[5 5 1]);
    hold off
    % pause(.00000001)
end

%% Initial condition
clear
clc

x_v = [1; 0; 0];
y_v = [0; 1; 0];
z_v = [0; 0; 1];
origin = [0; 0; 0];

PHI0 = 70;                                  % Initial roll angle (Between non-spinning frame and inertial) [deg]
PHI0 = deg2rad(PHI0);
PSI0 = 240;                                 % Initial yaw angle (Between non-spinning frame and inertial) [deg]
PSI0 = deg2rad(PSI0);
THETA0 = 0;                                 % Initial pitch angle (Between non-spinning frame and inertial) [deg]
THETA0 = deg2rad(THETA0);

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
T = RI01*RI02*RI03;
invTI0 = transpose(T);
phi0 = deg2rad(45);
theta0 = deg2rad(70);
psi0 = deg2rad(0);

R1 = [cos(phi0), sin(phi0), 0; -sin(phi0), cos(phi0), 0; 0, 0, 1];
R2 = [1, 0, 0; 0, cos(theta0), sin(theta0); 0, -sin(theta0), cos(theta0)];
R3 = [cos(psi0), sin(psi0), 0; -sin(psi0), cos(psi0), 0; 0, 0, 1];

% Combined rotation matrix
T = R3 * R2 * R1;   % Inertial to body
invT = transpose(T);% Body to inertial
invTI0 = T;
attx = invTI0*x_v;
atty = invTI0*y_v;
attz = invTI0*z_v;

figure(1)
plot3([origin(1) attx(1)],[origin(2) attx(2)],[origin(3) attx(3)],'r-^', 'LineWidth',3);
hold on
plot3([origin(1) atty(1)],[origin(2) atty(2)],[origin(3) atty(3)],'g-^', 'LineWidth',3);
hold on
plot3([origin(1) attz(1)],[origin(2) attz(2)],[origin(3) attz(3)],'b-^', 'LineWidth',3);
grid on;
xlabel('X axis'), ylabel('Y axis'), zlabel('Z axis')
legend('x_b', 'y_b', 'z_b')

%% 313 Rotation
clear
clc

x_v = [1; 0; 0];
y_v = [0; 1; 0];
z_v = [0; 0; 1];
origin = [0; 0; 0];
phi0 = 0;                                  % Initial roll angle (Between non-spinning frame and inertial) [deg]
phi0 = deg2rad(phi0);
psi0 = -30;                                 % Initial yaw angle (Between non-spinning frame and inertial) [deg]
psi0 = deg2rad(psi0);
theta0 = 70;                                 % Initial pitch angle (Between non-spinning frame and inertial) [deg]
theta0 = deg2rad(theta0);  

R1 = [cos(phi0), sin(phi0), 0; -sin(phi0), cos(phi0), 0; 0, 0, 1];
R2 = [1, 0, 0; 0, cos(theta0), sin(theta0); 0, -sin(theta0), cos(theta0)];
R3 = [cos(psi0), sin(psi0), 0; -sin(psi0), cos(psi0), 0; 0, 0, 1];

% Combined rotation matrix
T = R3 * R2 * R1;   % Inertial to body
invT = transpose(T);% Body to inertial

attx = T*x_v;
atty = T*y_v;
attz = T*z_v;

figure(1)
plot3([origin(1) attx(1)],[origin(2) attx(2)],[origin(3) attx(3)],'r-^', 'LineWidth',3);
hold on
plot3([origin(1) atty(1)],[origin(2) atty(2)],[origin(3) atty(3)],'g-^', 'LineWidth',3);
hold on
plot3([origin(1) attz(1)],[origin(2) attz(2)],[origin(3) attz(3)],'b-^', 'LineWidth',3);
grid on;
xlabel('X axis'), ylabel('Y axis'), zlabel('Z axis')
legend('x_b', 'y_b', 'z_b')

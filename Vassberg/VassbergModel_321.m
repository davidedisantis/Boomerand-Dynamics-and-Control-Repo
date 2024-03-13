%% Data
clear
close all
clc

tic
%--------------------------Simulation parameters--------------------------%
ts = 5;                     % Simulation time [s]

% Planetary parameters
g = 9.81;                % Gravity acceleration [m/s^2]
rho = 1.22;              % Air density at sea level [kg/m^3]
v_wind = [0; 0; 0];      % Wind speed in inertial reference frame [m/s]

%-----------------------Boomerang design parameters-----------------------%
l_blade = 3e-1; % Length of one blade [m]
nw = 2;         % Number of wings [ ]
R = 2.69e-1;    % Radius of rotation [m]
S = 2.28e-1;    % Disk area [m^2]
m = 1.30e-1;    % Boomerang mass [kg]
c = 4.88e-2;    % Mean chord [m]
S_b = l_blade*c;

x_ac = 7.45e-2; % Position of aerodynamic center in body coordinates
LAMBDAj = 120;  % Wing sweep angle [deg]
LAMBDAj = deg2rad(LAMBDAj);
gamma = 120;    % Folding angle [deg]
gamma = deg2rad(gamma);
betaj = 0;      % Wing coning angle 
betaj = deg2rad(betaj);
thetaj = 0;     % Wing pitch angle at root

%-----------------Moments of inertia of a single blade--------------------%
I_xi = 1.88e-3; 
I_eta = 4.78e-6;
I_zeta = 1.95e-3;
I_xieta = 5.14e-21;

Jj = [
    I_xi I_xieta 0;
    I_xieta I_eta 0;
    0 0 I_zeta
    ];

%-------------------------Launch conditions-------------------------------%
ThAng = 30;                                 % Throw angle from East direction
U0 = [25*cosd(ThAng); 25*sind(ThAng); 0];   % Initial throw speed in inertial frame [m/s]
omega0 = [0; 0; 10*pi*2];                   % Initial angular speed [rad/s]
R_pos0 = [0; 0; 1.5];                       % Initial position in inertial frame
PHI0 = 70;                           % Initial roll angle (Between non-spinning frame and inertial) [deg]
PHI0 = deg2rad(PHI0);
THETA0 = 0;                          % Initial pitch angle (Between non-spinning frame and inertial) [deg]
THETA0 = deg2rad(THETA0);
PSI0 = 240;                           % Initial yaw angle (Between non-spinning frame and inertial) [deg]
PSI0 = deg2rad(PSI0);
%------------Parameters for integration of aerodynamic forces-------------%
n = 50; % Number of intervals
l_integrate = linspace(0,l_blade,n+1);
h = l_blade/n;

%---------------------------Rotation matrices-----------------------------%
% From angular rates to Euler Angle rates

Tj = zeros(3, 3, nw);
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
    invTj(:,:,i+1) = transpose(Tj(:,:,i+1));
end
R1 = [1, 0, 0; 0, cos(PHI0), sin(PHI0); 0, -sin(PHI0), cos(PHI0)];
R2 = [cos(THETA0), 0, -sin(THETA0); 0, 1, 0; sin(THETA0), 0, cos(THETA0)];
R3 = [cos(PSI0), sin(PSI0), 0; -sin(PSI0), cos(PSI0), 0; 0, 0, 1];
% R1 = [cos(phi0), sin(phi0), 0; -sin(phi0), cos(phi0), 0; 0, 0, 1];
% R2 = [1, 0, 0; 0, cos(theta0), sin(theta0); 0, -sin(theta0), cos(theta0)];
% R3 = [cos(psi0), sin(psi0), 0; -sin(psi0), cos(psi0), 0; 0, 0, 1];

% Combined rotation matrix
% T0 = R3 * R2 * R1;   % Inertial to body
% invT0 = transpose(T0);% Body to inertial
T0 = R1*R2*R3;
invT0 = transpose(T0);

%------------Computation of moment of inertia wrt body frame--------------%
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

%-------------------Starting conditions in body frame---------------------%
% Eul0 = [phi0; theta0; psi0];
Eul0 = [PHI0; THETA0; PSI0];
u0 = T0*U0;
r0 = T0*R_pos0;
%------------------------Aerodynamic coefficients-------------------------%
Cl0 = 0.3484;
Cla = 5.7;
Cd0 = 0.02;
k = 1.5;
%-------------------ODE Function: definition and solving------------------%
f = @(t,s) [ ...
    ( FA(s(1:3), s(9), S_b, Cl0, Cla, Cd0, k, R, rho) + FG(s(4:6), m, g) )./m - cross(s(7:9),s(1:3)); ...
    EulRate321(s(4:5))*s(7:9);
    J\(MA(s(1:3), s(9), S_b, Cl0, Cla, R, rho) - cross(s(7:9), J*s(7:9)))];

[t,s] = ode45(f, [0 ts], [u0; Eul0; omega0]);

toc

%% Plots 
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

figure(1)
plot(t(cut),POS(:,3), '--k')
title('Altitude over time')
xlabel('Time [s]')
ylabel('Altitude [s]')

figure(2)
plot(POS(1,1), POS(1,2), 'go', ...
    POS(:,1), POS(:,2), 'k-', ...
    POS(end,1), POS(end,2), 'ro')
xlabel('Displacement along x [m]')
ylabel('Displacement along y [m]')
title('Trajectory - top view (xy plane)')
% hold on
% fimplicit(circle)
axis equal

figure(3)
plot(POS(1,1),POS(1,3), 'go', ...
    POS(:,1), POS(:,3), 'k-', ...
    POS(end,1), POS(end,3), 'ro')
xlabel('Displacement along x [m]')
ylabel('Displacement along z [m]')
title('Trajectory - xz plane')
axis equal

figure(4)
plot(POS(1,2), POS(1,3), 'go', ...
    POS(:,2), POS(:,3), 'k-', ...
    POS(end,2), POS(end,3), 'ro')
xlabel('Displacement along y [m]')
ylabel('Displacement along z [m]')
title('Trajectory - yz plane')
axis equal

figure(5)
plot(t(cut), vecnorm(u(cut,:)'))

% figure(6)
% plot3(POS(1,1), POS(1,2), POS(1,3), 'go', ...
%     POS(:,1), POS(:,2), POS(:,3),'k--', ...
%     POS(end,1), POS(end,2), POS(end,3), 'ro')
% grid on
% xlabel('X axis')
% ylabel('Y axis')
% zlabel('Z axis')
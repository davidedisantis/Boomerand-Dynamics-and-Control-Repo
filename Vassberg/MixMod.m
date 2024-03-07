clear
close all
clc

tic
% Simulation parameters
ts = .553;                     % Simulation time [s]

minstep = 1e-5;            % Min timestep for variable-step simulation [s]
maxstep = 1.66e-4;         % Max timestep for variable-step simulation [s]

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
S_b = R*c;

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
omega0 = [0; 0; 0*pi*2];            % Initial angular speed [rad/s]
PHI0 = 70;                           % Initial roll angle (Between non-spinning frame and inertial) [deg]
PHI0 = deg2rad(PHI0);
PSI0 = 240;                           % Initial yaw angle (Between non-spinning frame and inertial) [deg]
PSI0 = deg2rad(PSI0);
THETA0 = 0;                          % Initial pitch angle (Between non-spinning frame and inertial) [deg]
THETA0 = deg2rad(THETA0);  
R_pos0 = [0; 0; 1.5];                % Initial position in inertial frame

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

%------------------------Aerodynamic coefficients-------------------------%
Cl0 = 0.3484;
Cla = 5.2989;
%
f = @(t,s) [s(4:6); ...
    ( FA(s(4), s(5), s(6), s(12), S_b, Cl0, Cla, R, rho) + FG(s(7), s(8), s(9), m, g) )./m - cross(s(10:12),s(4:6)); ...
    s(10) + s(11)*sin(s(7))*tan(s(8)) + s(12)*cos(s(7))*tan(s(8)); ... 
    s(11)*cos(s(7)) - s(12)*sin(s(7)); ...
    s(11)*sin(s(7))/cos(s(8)) + s(12)*cos(s(7))/cos(s(8)); ...
    J\(MA(s(4), s(5), s(6), s(12), S_b, Cl0, Cla, R, rho) - cross(s(10:12), J*s(10:12)))];

[t,s] = ode45(f, [0 ts], [r0; u0; Eul0; omega0]);

toc

%%
x = s(:,1);
y = s(:,2);
z = s(:,3);
ux = s(:,4);
uy = s(:,5);
uz = s(:,6);

pos = [x, y,z];

phi = s(:,10);
theta = s(:,11);
psi = s(:,12);
POS = zeros(length(t), 3);
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

POS(i,:) = invTI0*invT0*pos(i,:)';
end
figure(1)
plot(t(1),POS(1,3), 'go', t,POS(:,3), '--k', t(end),POS(end,3), 'ro')

figure(2)
plot(POS(1,1), POS(1,2), 'go', ...
    POS(:,1), POS(:,2), 'k--', ...
    POS(end,1), POS(end,2), 'ro')
xlabel('Displacement along x [m]')
ylabel('Displacement along y [m]')
title('Trajectory - top view (xy plane)')

figure(3)
plot(POS(1,1),POS(1,3), 'go', ...
    POS(:,1), POS(:,3), 'k--', ...
    POS(end,1), POS(end,3), 'ro')
xlabel('Displacement along x [m]')
ylabel('Displacement along z [m]')
title('Trajectory - xz plane')

figure(4)
plot(POS(1,2), POS(1,3), 'go', ...
    POS(:,2), POS(:,3), 'k--', ...
    POS(end,2), POS(end,3), 'ro')
xlabel('Displacement along y [m]')
ylabel('Displacement along z [m]')
title('Trajectory - yz plane')

figure(5)
plot(t, sqrt(ux.^2 + uy.^2 + uz.^2))
function FG = FG(EulAng, m, g)
% R1 = [cos(phi), sin(phi), 0; -sin(phi), cos(phi), 0; 0, 0, 1];
% R2 = [1, 0, 0; 0, cos(theta), sin(theta); 0, -sin(theta), cos(theta)];
% R3 = [cos(psi), sin(psi), 0; -sin(psi), cos(psi), 0; 0, 0, 1];

phi = EulAng(1);
theta = EulAng(2);
psi = EulAng(3);

R01 = [
    1, 0, 0;
    0, cos(phi), sin(phi);
    0, -sin(phi), cos(phi)
    ];
R02 = [
    cos(theta), 0, -sin(theta);
    0, 1, 0;
    sin(theta), 0, cos(theta)
];
R03 = [
    cos(psi), sin(psi), 0;
    -sin(psi), cos(psi), 0;
    0, 0, 1
];
T = R01*R02*R03;

FG = -m*g*T*[0; 0; 1];
end
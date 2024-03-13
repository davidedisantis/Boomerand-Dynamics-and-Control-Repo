function A = EulRate321(ang)
phi = ang(1);
theta = ang(2);

A = [1 sin(phi)*tan(theta) cos(phi)*tan(theta); 
    0        cos(phi)         -sin(phi) 
    0    sin(phi)/cos(theta) cos(phi)/cos(theta)];
end
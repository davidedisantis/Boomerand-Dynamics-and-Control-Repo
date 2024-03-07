function A = EulRateMat(theta,psi)
A = [              sin(psi)/sin(theta), cos(psi)/sin(theta), 0;
                         cos(psi),          -sin(psi),       0;
-(cos(theta)*sin(psi))/sin(theta), -(cos(psi)*cos(theta))/sin(theta), 1];
end
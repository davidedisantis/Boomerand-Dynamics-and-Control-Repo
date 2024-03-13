function MA = MA(u, r, S_b, Cl0, Cla, R, rho)

chi = r*R/norm([u(1) u(2)]);

% alpha = atan2(-real(u(3)),norm([u(1) u(2)]));
alpha = -u(3)/norm([u(1) u(2)]);

% q = 0.5*rho*norm(u)^2;
q = 0.5*rho*norm([u(1) u(2)])^2;

M0 = 1/3*q*S_b*R*Cl0*chi;

if chi < 1
    Ma = (q*S_b*alpha*Cla*((pi-2*acos(-chi))/(16*pi*chi^2) + sqrt(1-chi^2)/(8*pi*chi) ...
    -(pi-2*acos(-chi))/(4*pi) + chi*sqrt(1-chi^2)/(4*pi)));
else 
    Ma = q*S_b*R*alpha*Cla*(1/4 - 1/(16*chi^2));
end
MA = 2*[M0 + Ma; 0; 0];

end
function MA = MA(u, r, S_b, Cl0, Cla, Cd0, Cm0, Cma, c, R, rho, nw)

chi = r*R/norm([u(1) u(2)]);

% alpha = atan2(-real(u(3)),norm([u(1) u(2)]));
alpha = -u(3)/norm([u(1) u(2)]);

% q = 0.5*rho*norm(u)^2;
q = 0.5*rho*norm([u(1) u(2)])^2;

Mx0 = 1/3*q*S_b*R*Cl0*chi;

if chi < 1
    Mxa = (q*S_b*alpha*Cla*((pi-2*acos(-chi))/(16*pi*chi^2) + sqrt(1-chi^2)/(8*pi*chi) ...
    -(pi-2*acos(-chi))/(4*pi) + chi*sqrt(1-chi^2)/(4*pi)));
else 
    Mxa = q*S_b*R*alpha*Cla*(1/4 - 1/(16*chi^2));
end

My0 = q*c*S_b*Cm0*(0.5 + 1/3*chi^2);
if chi < 1
    Mya = 1/pi*q*c*S_b*alpha*Cma*(3/2*sqrt(1-chi^2) + (chi + 1/(2*chi))*(acos(-chi)-pi/2));
else
    Mya = q*c*S_b*alpha*Cma*(chi/2 + 1/(4*chi));
end

Mx = Mx0 + Mxa;
My = My0 + Mya;
Mz = 1/3*q*S_b*R*Cd0*chi;

MA = nw*[Mx; My; Mz];

end
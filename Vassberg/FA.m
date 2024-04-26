function FA = FA(u, r, S_b, Cl0, Cla, Cd0, k, R, rho, nw)

chi = r*R/norm([u(1) u(2)]);

alpha = -u(3)/norm([u(1) u(2)]);

q = 0.5*rho*norm(u)^2;

L0 = q*S_b*Cl0*(0.5 + 1/3*chi^2);
if chi < 1
    La = 1/pi*q*S_b*alpha*Cla*(3/2*sqrt(1-chi^2) + (chi + 1/(2*chi))*(acos(-chi)-pi/2));
else
    La = q*S_b*alpha*Cla*(chi/2 + 1/(4*chi));
end
L = L0 + La;

D0 = q*S_b*Cd0*(1/2 + 1/3*chi^2);
Da = q*S_b*k*alpha^2;
D = D0 + Da;

beta = atan2(u(2),u(1));
FA = nw*[-D*cos(beta); - D*sin(beta); L];

end
function FA = FA(u, r, S_b, Cl0, Cla, R, rho)

chi = r*R/norm(u);

alpha = atan2(-u(3),norm([u(1) u(2)]));

q = 0.5*rho*norm(u)^2;

L0 = q*S_b*Cl0*(0.5 + 1/3*chi^2);
if chi < 1
    La = 1/pi*q*S_b*alpha*Cla*(3/2*sqrt(1-chi^2) + (chi + 1/(2*chi))*(acos(-chi)-pi/2));
else
    La = q*S_b*alpha*Cla*(chi/2 + 1/(4*chi));
end

beta = atan2(real(u(2)),real(u(1)));

FA = 2*[-(L0+La)*sin(alpha)*cos(beta); -(L0+La)*sin(alpha)*sin(beta); (L0 + La)];

end
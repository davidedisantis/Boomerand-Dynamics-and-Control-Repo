function y = myMeasurementFcn(x)
% States: [pos(3x1), vel(3x1), acc(3x1), EulAng(3x1), omega(3x1),
% omegadot(3x1)
% IMU only measures acc and omega

y = [x(7:9) x(13:15)];
end
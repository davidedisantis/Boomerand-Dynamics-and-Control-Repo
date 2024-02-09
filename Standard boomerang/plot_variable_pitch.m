%% Plot results with pitch variation markers - Use this section if no comparison with constant pitch is needed

R_inertial = zeros(length(sim.tout),3);
u = zeros(length(sim.tout),3);
normU = zeros(length(sim.tout),1);

for i = 1:length(sim.tout)
    R_inertial(i, :) = sim.R_inertial(:,1,i);
    u(i,:) = sim.U(:,1,i);
    normU(i) = norm(u(i,:));
end

act_ind = ceil(length(sim.tout)*t_action/t); % Index of time instant at which pitch variation (action) happens

% Plots
figure(1)
plot(sim.tout, R_inertial(:,3), 'k-')
xlabel('Time [s]')
ylabel('Altitude [m]')

figure(2)
plot(R_inertial(1,1), R_inertial(1,2), 'go', ...
    R_inertial(:,1), R_inertial(:,2), 'k-', ...
    R_inertial(act_ind,1), R_inertial(act_ind,2), 'bo', ...
    R_inertial(end,1), R_inertial(end,2), 'ro')
xlabel('Displacement along x [m]')
ylabel('Displacement along y [m]')
title('Trajectory - top view (xy plane)')
legend('Starting point','Trajectory', 'Change of pitch', 'End point')

figure(3)
plot(R_inertial(1,1), R_inertial(1,3), 'go', ...
    R_inertial(:,1), R_inertial(:,3), 'k-', ...
    R_inertial(act_ind,1), R_inertial(act_ind,3), 'bo', ...
    R_inertial(end,1), R_inertial(end,3), 'ro')
xlabel('Displacement along x [m]')
ylabel('Displacement along z [m]')
title('Trajectory - xz plane')
legend('Starting point','Trajectory', 'Change of pitch', 'End point')

figure(4)
plot(R_inertial(1,2), R_inertial(1,3), 'go', ...
    R_inertial(:,2), R_inertial(:,3), 'k-', ...
    R_inertial(act_ind,2), R_inertial(act_ind,3), 'bo', ...
    R_inertial(end,2), R_inertial(end,3), 'ro')
xlabel('Displacement along y [m]')
ylabel('Displacement along z [m]')
title('Trajectory - yz plane')
legend('Starting point','Trajectory', 'Change of pitch', 'End point')

figure(5)
plot(sim.tout, normU, 'k-')
xlabel('Time [s]')
ylabel('Speed U [m/s]')

figure(6)
plot3(R_inertial(1,1), R_inertial(1,2), R_inertial(1,3), 'go', ...
    R_inertial(:,1), R_inertial(:,2), R_inertial(:,3), 'k-', ...
    R_inertial(act_ind,1), R_inertial(act_ind,2), R_inertial(act_ind,3), 'bo', ...
    R_inertial(end,1), R_inertial(end,2), R_inertial(end,3), 'ro')
title('3D Trajectory')
legend('Starting point','Trajectory', 'Change of pitch', 'End point')
grid on

%% Plot results with pitch variation markers - No legends

R_inertial = zeros(length(sim.tout),3);
u = zeros(length(sim.tout),3);
normU = zeros(length(sim.tout),1);

for i = 1:length(sim.tout)
    R_inertial(i, :) = sim.R_inertial(:,1,i);
    u(i,:) = sim.U(:,1,i);
    normU(i) = norm(u(i,:));
end

act_ind = ceil(length(sim.tout)*t_action/t); % Index of time instant at which pitch variation (action) happens

% Plots
figure(1);
hold on
plot(sim.tout, R_inertial(:,3), 'k-')
xlabel('Time [s]')
ylabel('Altitude [m]')

figure(2);
hold on
plot(R_inertial(1,1), R_inertial(1,2), 'go', ...
    R_inertial(:,1), R_inertial(:,2), 'k-', ...
    R_inertial(act_ind,1), R_inertial(act_ind,2), 'bo', ...
    R_inertial(end,1), R_inertial(end,2), 'ro')
xlabel('Displacement along x [m]')
ylabel('Displacement along y [m]')
title('Trajectory - top view (xy plane)')

figure(3);
hold on
plot(R_inertial(1,1), R_inertial(1,3), 'go', ...
    R_inertial(:,1), R_inertial(:,3), 'k-', ...
    R_inertial(act_ind,1), R_inertial(act_ind,3), 'bo', ...
    R_inertial(end,1), R_inertial(end,3), 'ro')
xlabel('Displacement along x [m]')
ylabel('Displacement along z [m]')
title('Trajectory - xz plane')

figure(4);
hold on
plot(R_inertial(1,2), R_inertial(1,3), 'go', ...
    R_inertial(:,2), R_inertial(:,3), 'k-', ...
    R_inertial(act_ind,2), R_inertial(act_ind,3), 'bo', ...
    R_inertial(end,2), R_inertial(end,3), 'ro')
xlabel('Displacement along y [m]')
ylabel('Displacement along z [m]')
title('Trajectory - yz plane')

figure(5);
hold on
plot(sim.tout, normU, 'k-')
xlabel('Time [s]')
ylabel('Speed U [m/s]')

figure(6);
hold on
plot3(R_inertial(1,1), R_inertial(1,2), R_inertial(1,3), 'go', ...
    R_inertial(:,1), R_inertial(:,2), R_inertial(:,3), 'k-', ...
    R_inertial(act_ind,1), R_inertial(act_ind,2), R_inertial(act_ind,3), 'bo', ...
    R_inertial(end,1), R_inertial(end,2), R_inertial(end,3), 'ro')
title('3D Trajectory')
grid on

%% Plot results with pitch variation markers - Run only after plotting no pitch variation results
% Run the code, plot the results by running prev section, comment "close all",
% re-run the code, run this section
R_inertial = zeros(length(sim.tout),3);
u = zeros(length(sim.tout),3);
normU = zeros(length(sim.tout),1);

for i = 1:length(sim.tout)
    R_inertial(i, :) = sim.R_inertial(:,1,i);
    u(i,:) = sim.U(:,1,i);
    normU(i) = norm(u(i,:));
end

act_ind = ceil(length(sim.tout)*t_action/t); % Index of time instant at which pitch variation (action) happens

% Plots
figure(1);
hold on
plot(sim.tout, R_inertial(:,3), 'k--')
xlabel('Time [s]')
ylabel('Altitude [m]')
legend('Constant pitch', 'Variable pitch')

figure(2);
hold on
plot(R_inertial(:,1), R_inertial(:,2), 'k--', ...
    R_inertial(act_ind,1), R_inertial(act_ind,2), 'bo', ...
    R_inertial(end,1), R_inertial(end,2), 'ro')
xlabel('Displacement along x [m]')
ylabel('Displacement along y [m]')
title('Trajectory - top view (xy plane)')
legend('Starting point', 'Constant pitch', '', 'Variable pitch', 'Change of pitch', 'End point')

figure(3);
hold on
plot(R_inertial(:,1), R_inertial(:,3), 'k--', ...
    R_inertial(act_ind,1), R_inertial(act_ind,3), 'bo', ...
    R_inertial(end,1), R_inertial(end,3), 'ro')
xlabel('Displacement along x [m]')
ylabel('Displacement along z [m]')
title('Trajectory - xz plane')
legend('Starting point', 'Constant pitch', '', 'Variable pitch', 'Change of pitch', 'End point')

figure(4);
hold on
plot(R_inertial(:,2), R_inertial(:,3), 'k--', ...
    R_inertial(act_ind,2), R_inertial(act_ind,3), 'bo', ...
    R_inertial(end,2), R_inertial(end,3), 'ro')
xlabel('Displacement along y [m]')
ylabel('Displacement along z [m]')
title('Trajectory - yz plane')
legend('Starting point', 'Constant pitch', '', 'Variable pitch', 'Change of pitch', 'End point')

figure(5);
hold on
plot(sim.tout, normU, 'k--')
xlabel('Time [s]')
ylabel('Speed U [m/s]')
legend('Constant pitch', 'Variable pitch')

figure(6);
hold on
plot3(R_inertial(:,1), R_inertial(:,2), R_inertial(:,3), 'k--', ...
    R_inertial(act_ind,1), R_inertial(act_ind,2), R_inertial(act_ind,3), 'bo', ...
    R_inertial(end,1), R_inertial(end,2), R_inertial(end,3), 'ro')
title('3D Trajectory')
legend('Starting point', 'Constant pitch', '', 'Variable pitch', 'Change of pitch', 'End point')
%% Plot results over variable pitch plots

R_inertial = zeros(length(sim.tout),3);
u = zeros(length(sim.tout),3);
normU = zeros(length(sim.tout),1);

for i = 1:length(sim.tout)
    R_inertial(i, :) = sim.R_inertial(:,1,i);
    u(i,:) = sim.U(:,1,i);
    normU(i) = norm(u(i,:));
end
% Identify index corresponding to desired time instants - allows to find
% the indexes relative to a desired time instant even for variable step sim
[~, act_ind] = min(abs(sim.tout - t_action));
[~, act_stop] = min(abs(sim.tout) - (t_action+t_rise));

% Plots
figure(1);
hold on
plot(sim.tout, R_inertial(:,3), 'k--')
xlabel('Time [s]')
ylabel('Altitude [m]')
legend('Step-variation in pitch', 'Gradual variation')
title('Altitude over time')

figure(2);
hold on
plot(R_inertial(:,1), R_inertial(:,2), 'k--', ...
    R_inertial(act_ind,1), R_inertial(act_ind,2), 'bo', ...
    R_inertial(act_stop,1), R_inertial(act_stop,2), 'co', ...
    R_inertial(end,1), R_inertial(end,2), 'ro')
xlabel('Displacement along x [m]')
ylabel('Displacement along y [m]')
title('Trajectory - top view (xy plane)')
legend('Starting point', 'Step-variation in pitch', '', '', 'Gradual pitch variation')

figure(3);
hold on
plot(R_inertial(1,1), R_inertial(1,3), 'go', ...
    R_inertial(:,1), R_inertial(:,3), 'k--', ...
    R_inertial(act_ind,1), R_inertial(act_ind,3), 'bo', ...
    R_inertial(act_stop,1), R_inertial(act_stop,3), 'co', ...
    R_inertial(end,1), R_inertial(end,3), 'ro')
xlabel('Displacement along x [m]')
ylabel('Displacement along z [m]')
title('Trajectory - xz plane')

figure(4);
hold on
plot(R_inertial(1,2), R_inertial(1,3), 'go', ...
    R_inertial(:,2), R_inertial(:,3), 'k--', ...
    R_inertial(act_ind,2), R_inertial(act_ind,3), 'bo', ...
    R_inertial(act_stop,2), R_inertial(act_stop,3), 'co', ...
    R_inertial(end,2), R_inertial(end,3), 'ro')
xlabel('Displacement along y [m]')
ylabel('Displacement along z [m]')
title('Trajectory - yz plane')

figure(5);
hold on
plot(sim.tout, normU, 'k--')
xlabel('Time [s]')
ylabel('Speed U [m/s]')

figure(6);
hold on
plot3(R_inertial(1,1), R_inertial(1,2), R_inertial(1,3), 'go', ...
    R_inertial(:,1), R_inertial(:,2), R_inertial(:,3), 'k--', ...
    R_inertial(act_ind,1), R_inertial(act_ind,2), R_inertial(act_ind,3), 'bo', ...
    R_inertial(act_stop,1), R_inertial(act_stop,2), R_inertial(act_stop,3), 'co', ...
    R_inertial(end,1), R_inertial(end,2), R_inertial(end,3), 'ro')
title('3D Trajectory')
grid on
%% Plot results with pitch variation markers - No legends - Two actuations - First print

R_inertial = zeros(length(sim.tout),3);
u = zeros(length(sim.tout),3);
normU = zeros(length(sim.tout),1);

for i = 1:length(sim.tout)
    R_inertial(i, :) = sim.R_inertial(:,1,i);
    u(i,:) = sim.U(:,1,i);
    normU(i) = norm(u(i,:));
end
% Identify index corresponding to desired time instants - allows to find
% the indexes relative to a desired time instant even for variable step sim
[~, act_ind1] = min(abs(sim.tout - t_action1));
[~, act_ind2] = min(abs(sim.tout - t_action2));

% Plots
figure(1)
plot(sim.tout, R_inertial(:,3), 'k-')
xlabel('Time [s]')
ylabel('Altitude [m]')

figure(2)
plot(R_inertial(1,1), R_inertial(1,2), 'go', ...
    R_inertial(:,1), R_inertial(:,2), 'k-', ...
    R_inertial(act_ind1,1), R_inertial(act_ind1,2), 'bo', ...
    R_inertial(act_ind2,1), R_inertial(act_ind2,2), 'co', ...
    R_inertial(end,1), R_inertial(end,2), 'ro')
xlabel('Displacement along x [m]')
ylabel('Displacement along y [m]')
title('Trajectory - top view (xy plane)')

figure(3)
plot(R_inertial(1,1), R_inertial(1,3), 'go', ...
    R_inertial(:,1), R_inertial(:,3), 'k-', ...
    R_inertial(act_ind1,1), R_inertial(act_ind1,3), 'bo', ...
    R_inertial(act_ind2,1), R_inertial(act_ind2,3), 'co', ...
    R_inertial(end,1), R_inertial(end,3), 'ro')
xlabel('Displacement along x [m]')
ylabel('Displacement along z [m]')
title('Trajectory - xz plane')

figure(4)
plot(R_inertial(1,2), R_inertial(1,3), 'go', ...
    R_inertial(:,2), R_inertial(:,3), 'k-', ...
    R_inertial(act_ind1,2), R_inertial(act_ind1,3), 'bo', ...
    R_inertial(act_ind2,2), R_inertial(act_ind2,3), 'co', ...
    R_inertial(end,2), R_inertial(end,3), 'ro')
xlabel('Displacement along y [m]')
ylabel('Displacement along z [m]')
title('Trajectory - yz plane')

figure(5)
plot(sim.tout, normU, 'k-')
xlabel('Time [s]')
ylabel('Speed U [m/s]')

figure(6)
plot3(R_inertial(1,1), R_inertial(1,2), R_inertial(1,3), 'go', ...
    R_inertial(:,1), R_inertial(:,2), R_inertial(:,3), 'k-', ...
    R_inertial(act_ind1,1), R_inertial(act_ind1,2), R_inertial(act_ind1,3), 'bo', ...
    R_inertial(act_ind2,1), R_inertial(act_ind2,2), R_inertial(act_ind2,3), 'co', ...
    R_inertial(end,1), R_inertial(end,2), R_inertial(end,3), 'ro')
title('3D Trajectory')
grid on
%% Plot results with pitch variation markers - No legends - Two actuations - Second print

R_inertial = zeros(length(sim.tout),3);
u = zeros(length(sim.tout),3);
normU = zeros(length(sim.tout),1);

for i = 1:length(sim.tout)
    R_inertial(i, :) = sim.R_inertial(:,1,i);
    u(i,:) = sim.U(:,1,i);
    normU(i) = norm(u(i,:));
end
% Identify index corresponding to desired time instants - allows to find
% the indexes relative to a desired time instant even for variable step sim
[~, act_ind1] = min(abs(sim.tout - t_action1));
[~, act_ind2] = min(abs(sim.tout - t_action2));

% Plots
figure(1);
hold on
plot(sim.tout, R_inertial(:,3), 'k-')
xlabel('Time [s]')
ylabel('Altitude [m]')

figure(2);
hold on
plot(R_inertial(1,1), R_inertial(1,2), 'go', ...
    R_inertial(:,1), R_inertial(:,2), 'k-', ...
    R_inertial(act_ind1,1), R_inertial(act_ind1,2), 'bo', ...
    R_inertial(act_ind2,1), R_inertial(act_ind2,2), 'co', ...
    R_inertial(end,1), R_inertial(end,2), 'ro')
xlabel('Displacement along x [m]')
ylabel('Displacement along y [m]')
title('Trajectory - top view (xy plane)')

figure(3);
hold on
plot(R_inertial(1,1), R_inertial(1,3), 'go', ...
    R_inertial(:,1), R_inertial(:,3), 'k-', ...
    R_inertial(act_ind1,1), R_inertial(act_ind1,3), 'bo', ...
    R_inertial(act_ind2,1), R_inertial(act_ind2,3), 'co', ...
    R_inertial(end,1), R_inertial(end,3), 'ro')
xlabel('Displacement along x [m]')
ylabel('Displacement along z [m]')
title('Trajectory - xz plane')

figure(4);
hold on
plot(R_inertial(1,2), R_inertial(1,3), 'go', ...
    R_inertial(:,2), R_inertial(:,3), 'k-', ...
    R_inertial(act_ind1,2), R_inertial(act_ind1,3), 'bo', ...
    R_inertial(act_ind2,2), R_inertial(act_ind2,3), 'co', ...
    R_inertial(end,2), R_inertial(end,3), 'ro')
xlabel('Displacement along y [m]')
ylabel('Displacement along z [m]')
title('Trajectory - yz plane')

figure(5);
hold on
plot(sim.tout, normU, 'k-')
xlabel('Time [s]')
ylabel('Speed U [m/s]')

figure(6);
hold on
plot3(R_inertial(1,1), R_inertial(1,2), R_inertial(1,3), 'go', ...
    R_inertial(:,1), R_inertial(:,2), R_inertial(:,3), 'k-', ...
    R_inertial(act_ind1,1), R_inertial(act_ind1,2), R_inertial(act_ind1,3), 'bo', ...
    R_inertial(act_ind2,1), R_inertial(act_ind2,2), R_inertial(act_ind2,3), 'co', ...
    R_inertial(end,1), R_inertial(end,2), R_inertial(end,3), 'ro')
title('3D Trajectory')
grid on
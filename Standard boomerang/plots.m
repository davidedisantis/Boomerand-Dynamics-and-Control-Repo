u = sim.U(:,:);
normU = vecnorm(u);
PHI = sim.PHI;
load('PHI_natural')

figure(1)
plot(sim.tout, normU, 'k-')
xlabel('Time [s]')
ylabel('Speed U [m/s]')
% xlim([0 .3])
figure(2)
plot(sim.tout,rad2deg(PHI_des)*ones(size(sim.tout)), sim.tout, rad2deg(PHI), 'LineWidth',1);
hold on
plot(t1, PHI1, 'k--', 'LineWidth', 0.5)
xlim([0, sim.tout(end)])
legend('Desired angle', 'Actual angle', 'Uncontrolled angle - natural trajectory', 'FontSize', 11)
% xlim([0 .3])

% Averages
% n = floor(length(sim.tout)/500);
% n1 = floor(length(t1) / 500);
% % Truncate the vector to contain only full groups of 10 elements
% truncated_PHI = PHI(1:n * 500);
% truncated_t = sim.tout(1:n*500);
% truncated_PHI1 = PHI1(1:n1 * 500);
% truncated_t1 = t1(1:n1*500);
% % Reshape the truncated vector into a matrix
% matrixt = reshape(truncated_t, 500, n);
% matrixPHI = reshape(truncated_PHI, 500, n);
% matrixt1 = reshape(truncated_t1, 500, n1);
% matrixPHI1 = reshape(truncated_PHI1, 500, n1);
% % Calculate the mean along one dimension (e.g., rows)
% mean_t = mean(matrixt, 1); % '1' because of column-wise means
% mean_PHI = mean(matrixPHI, 1); % '1' because of column-wise means
% mean_t1 = mean(matrixt1, 1); % '1' because of column-wise means
% mean_PHI1 = mean(matrixPHI1, 1); % '1' because of column-wise means
figure(3)
plot(sim.tout,rad2deg(PHI_des)*ones(size(sim.tout)), mean_t, rad2deg(mean_PHI), 'LineWidth',1);
hold on
plot(mean_t1, mean_PHI1, 'k--','LineWidth', 0.5)
xlim([0, sim.tout(end)])
legend('Desired angle', 'Avg Actual angle', 'Avg Uncontrolled angle - natural trajectory', 'FontSize', 11)

%% Trajectory and attitude of the non-spinning disk
phi = sim.EulAng(:,1);
theta = sim.EulAng(:,2);
psi = sim.EulAng(:,3);
PHI = sim.EULANG(1,:);
THETA = sim.EULANG(2,:);
PSI = sim.EULANG(3,:);
R_inertial = sim.R_inertial(:,:)';
% % Define the filename for the video
% videoFilename = 'trajectory_video.avi';

% % Create a VideoWriter object
% video = VideoWriter(videoFilename);
% 
% % Set the frame rate (frames per second)
% numFrames = length(phi(1:50:end));
% frameRate = numFrames/(sim.tout(end)); % Adjust as needed
% video.FrameRate = frameRate;
% 
% % Open the VideoWriter
% open(video);

x_v = [1; 0; 0];
y_v = [0; 1; 0];
z_v = [0; 0; 1];

for i = 1:50:length(sim.tout)
    R1 = [1, 0, 0; 0, cos(phi(i)), sin(phi(i)); 0, -sin(phi(i)), cos(phi(i))];
    R2 = [cos(theta(i)), 0, -sin(theta(i)); 0, 1, 0; sin(theta(i)), 0, cos(theta(i))];
    R3 = [cos(psi(i)), sin(psi(i)), 0; -sin(psi(i)), cos(psi(i)), 0; 0, 0, 1];
    T0 = R1*R2*R3;
    invT0 = transpose(T0);

    attx = invTI0*invT0*x_v;
    atty = invTI0*invT0*y_v;
    attz = invTI0*invT0*z_v;

    R1n = [1, 0, 0; 0, cos(PHI(i)), sin(PHI(i)); 0, -sin(PHI(i)), cos(PHI(i))];
    R2n = [cos(THETA(i)), 0, -sin(THETA(i)); 0, 1, 0; sin(THETA(i)), 0, cos(THETA(i))];
    R3n = [cos(PSI(i)), sin(PSI(i)), 0;-sin(PSI(i)), cos(PSI(i)), 0; 0, 0, 1];
    T = R1n*R2n*R3n;

    invT = transpose(T);
    attxn = invTI0*invT*x_v;
    attyn = invTI0*invT*y_v;
    attzn = invTI0*invT*z_v;

    figure(1);
    plot3(R_inertial(:,1),R_inertial(:,2),R_inertial(:,3),'k-')
    hold on
    quiver3(R_inertial(i,1),R_inertial(i,2),R_inertial(i,3), attx(1), attx(2), attx(3), 'b')
    hold on
    quiver3(R_inertial(i,1),R_inertial(i,2),R_inertial(i,3), atty(1), atty(2), atty(3), 'b')
    hold on
    quiver3(R_inertial(i,1),R_inertial(i,2),R_inertial(i,3), attz(1), attz(2), attz(3), 'r')
    hold on
    quiver3(R_inertial(i,1),R_inertial(i,2),R_inertial(i,3), attxn(1), attxn(2), attxn(3), 'c')
    hold on
    quiver3(R_inertial(i,1),R_inertial(i,2),R_inertial(i,3), attyn(1), attyn(2), attyn(3), 'c')
    hold on
    quiver3(R_inertial(i,1),R_inertial(i,2),R_inertial(i,3), attzn(1), attzn(2), attzn(3), 'm')    
    xlabel('X axis'), ylabel('Y axis'), zlabel('Z axis')
    hold off
    legend('Trajectory', 'x_b', 'y_b', 'z_b','interpreter', 'TeX')
    grid on

    xlim([min(R_inertial(:,1))-1 max(R_inertial(:,1))+1])
    ylim([min(R_inertial(:,2))-1 max(R_inertial(:,2))+1])
    zlim([min(R_inertial(:,3))-1 max(R_inertial(:,3))+1])
    %  % Capture the current frame
    % frame = getframe(gcf);
    % 
    % % Write the frame to the video
    % writeVideo(video, frame);
end
% % Close the VideoWriter
% close(video);
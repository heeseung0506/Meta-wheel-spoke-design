% Force–Local_Deflection_Curve_for_Inclined_Cantilever_Beam

% Beam properties
E = 7.4e6;      % Young's modulus (Pa)
I = 1e-6;       % Second moment of area (m^4)
L = 1.0;        % Beam length (m)

% Angles in degrees
theta_list = [0, 10, 20, 30];
theta_rad_list = deg2rad(theta_list);

% Force range (N)
P = linspace(0, 500, 300);  % Applied vertical force

% Initialize figure
figure; hold on;

% Loop over each angle
for i = 1:length(theta_list)
    theta = theta_rad_list(i);
    P_perp = P .* cos(theta);                          % Perpendicular force
    delta_local = (P_perp * L^3) / (3 * E * I);        % Deflection (m)
    delta_local_mm = delta_local * 1e3;                % Convert to mm
    plot(delta_local_mm, P, 'LineWidth', 2, ...
         'DisplayName', sprintf('\\beta = %d^\\circ', theta_list(i)));
end

% Plot settings
xlabel('Local Deflection \delta_{local} (mm)', 'FontSize', 12);
ylabel('Applied Load P(N)', 'FontSize', 12);
title('Force–Local Deflection Curve for Inclined Cantilever Beam', 'FontSize', 13);
legend('Location', 'northwest');
grid on;
xlim([0 inf]);

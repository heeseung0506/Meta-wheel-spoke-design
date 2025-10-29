% Inclined Cantilever Beam Deflection Visualization

% Parameters
P = 100;         % Applied load (N)
L = 1.0;         % Beam length (m)
E = 200e9;       % Young's modulus (Pa)
I = 1e-6;        % Second moment of area (m^4)

% Angle range (degrees)
theta_deg = linspace(0, 45, 300);     % from 0° to 90°
theta_rad = deg2rad(theta_deg);       % convert to radians

% Deflections
v_horizontal = P*L^3 / (3*E*I);                         % baseline (horizontal beam)
delta_local = (P .* cos(theta_rad)) .* L^3 ./ (3*E*I);  % local deflection
v_global = delta_local .* cos(theta_rad);               % projected onto global vertical

% Plot
figure;
plot(theta_deg, v_global*1e3, 'b--', 'LineWidth', 2); hold on;
plot(theta_deg, delta_local*1e3, 'r--', 'LineWidth', 2);
yline(v_horizontal*1e3, 'g-', 'LineWidth', 3);

xlabel('\theta (degrees)');
ylabel('Deflection (mm)');
title('Deflection of Inclined Cantilever Beam under Vertical Load');
legend('Global Vertical Deflection v_{max}^{(\theta)}', ...
       'Local Deflection \delta_{local}', ...
       'Cantilever Beam Reference', ...
       'Location', 'southwest');

grid on;

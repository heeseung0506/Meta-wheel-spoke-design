% Load data from CSV files
data1 = readmatrix('contin.csv'); 
data2 = readmatrix('0degree.csv'); 
data3 = readmatrix('10degree.csv'); 
data4 = readmatrix('20degree.csv'); 
data5 = readmatrix('2DTWEEL2.csv');

% Extract angles and forces
angles1 = data1(:,1); % First column: Angles in radians
adjusted_force1 = data1(:,2); % Second column: Adjusted force values

angles2 = data2(:,1);
adjusted_force2 = data2(:,2);

angles3 = data3(:,1);
adjusted_force3 = data3(:,2);

angles4 = data4(:,1);
adjusted_force4 = data4(:,2);

angles5 = data5(:,1);
adjusted_force5 = data5(:,2);

% Rotate angles by 180 degrees
angles_rotated1 = mod(angles1 + pi, 2*pi);
angles_rotated2 = mod(angles2 + pi, 2*pi);
angles_rotated3 = mod(angles3 + pi, 2*pi);
angles_rotated4 = mod(angles4 + pi, 2*pi);
angles_rotated5 = mod(angles5 + pi, 2*pi);

% Create polar plot
figure;

% Plot each dataset with different markers
polarplot(angles_rotated1, adjusted_force1, '-go', 'Marker', 'o', 'LineWidth', 1.2); 
hold on;
polarplot(angles_rotated2, adjusted_force2, '-rs', 'Marker', 's', 'LineWidth', 1.2); % 네모
polarplot(angles_rotated3, adjusted_force3, '-b^', 'Marker', '^', 'LineWidth', 1.2); % 세모
polarplot(angles_rotated4, adjusted_force4, '-kd', 'Marker', 'd', 'LineWidth', 1.2); % 마름모
polarplot(angles_rotated5, adjusted_force5, '-cx', 'Marker', 'x', 'LineWidth', 1.2); 


% Customize the plot
ax = gca;
ax.ThetaZeroLocation = 'top'; % Set the zero location at the top
ax.ThetaDir = 'clockwise'; % Set direction to clockwise

% Add title and legend
title('Radial Force Distribution','FontSize', 16);
legend('Continuous', '0 degree', '10 degree', '20 degree', '2D TWEEL','FontSize', 14);

% Show the plot
hold off;

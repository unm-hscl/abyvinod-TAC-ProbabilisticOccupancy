% Define the permitted intercept zones for Figure 4
clf

% Define polytope for area of interest
area_of_interest = Polyhedron('V', ...
    [0, 0;              % University and Central
     5, -0.125;         % Terrace and Central
     35,-0.025;          % Central and Girard
     35.2, 10;           % Campus and Girard
     20.5, 13.5;         % Las Lomas and Campus
     1, 13]);            % Las Lomas and University

% Define polytopes for the left fly zones
pursuer_cvx = [ones(2,0)*Polyhedron()];
pursuer_cvx(1,1) = Polyhedron('V', ...
    [0, 0;              % University and Central
     1, 13;              % Las Lomas and University
     2.6, 13.04;         % Redondo and Las Lomas
     1.7, 1]);        % Bend of Redondo
pursuer_cvx(1,2) = Polyhedron('V', ...
    [0, 0;              % University and Central
     1.7,  1;           % Bend of Redondo
     15,  0.85;         % Harvard into the University
     15,  -0.0925;      % Harvard and Central
     5, -0.125]);       % Terrace and Central
pursuer_cvx(1,3) = ones(2,0) * Polyhedron();

% Define polytopes for the top fly zones
pursuer_cvx(2,1) = Polyhedron('V', ...
    [8, 13.18;
     15, 13.358;
     15.1, 9.5;
     8, 9]);
pursuer_cvx(2,2) = Polyhedron('V', ...
    [8, 9;
     7.25,7.8;
     13.65, 8.75;
     13.5, 9.75]);
pursuer_cvx(2,3) = Polyhedron('V', ...
    [13.25, 9;
     14.5, 8;
     18.5, 8;
     18, 9]);

% Define polytopes for the bottom fly zones
pursuer_cvx(3,1) = Polyhedron('V', ...
    [35.025, 1.25;
     35,-0.025;          % Central and Girard
     22.5, -0.0655;
     22.5, 1.25]);
pursuer_cvx(3,2) = Polyhedron('V', ...
    [25, 1.25;
     30, 1.25;
     30, 7.3;
     25.1, 7.3]);
pursuer_cvx(3,3) = Polyhedron('V', ...
    [25.5, 6.15;
     25.5, 9;
     20.5, 6.15;
     21.5, 9]);
 
% Plot
plot(area_of_interest, 'color', 'b', 'alpha',0.2);
hold on;
plot(pursuer_cvx(1,1), 'color', 'g', 'alpha',0.4);
plot(pursuer_cvx(1,2), 'color', 'g', 'alpha',0.4);
plot(pursuer_cvx(1,3), 'color', 'g', 'alpha',0.4);
plot(pursuer_cvx(2,1), 'color', 'g', 'alpha',0.4);
plot(pursuer_cvx(2,2), 'color', 'g', 'alpha',0.4);
plot(pursuer_cvx(2,3), 'color', 'g', 'alpha',0.4);
plot(pursuer_cvx(3,1), 'color', 'g', 'alpha',0.4);
plot(pursuer_cvx(3,2), 'color', 'g', 'alpha',0.4);
plot(pursuer_cvx(3,3), 'color', 'g', 'alpha',0.4);
axis equal;
box on;
set(gca,'FontSize',fontSize*2);
grid on;
xlim([-2,37]);
ylim([-2,15]);
xlabel('x','Interpreter','latex');
ylabel('y','Interpreter','latex');
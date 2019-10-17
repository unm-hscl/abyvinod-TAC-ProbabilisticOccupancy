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
pursuer1_cvx1 = Polyhedron('V', ...
    [0, 0;              % University and Central
     1, 13;              % Las Lomas and University
     2.6, 13.04;         % Redondo and Las Lomas
     1.7, 0.5]);        % Bend of Redondo
pursuer1_cvx2 = Polyhedron('V', ...
    [0, 0;              % University and Central
     1.7,  0.5;         % Bend of Redondo
     15,  0.25;         % Harvard into the University
     15,  -0.0925;      % Harvard and Central
     5, -0.125]);       % Terrace and Central

% Define polytopes for the top fly zones
pursuer2_cvx1 = Polyhedron('V', ...
    [8, 13.18;
     15, 13.358;
     15.1, 9.5;
     8, 9]);
pursuer2_cvx2 = Polyhedron('V', ...
    [8, 9;
     7.25,7.8;
     13.65, 8.75;
     13.5, 9.75]);
pursuer2_cvx3 = Polyhedron('V', ...
    [13.25, 9;
     14.5, 8;
     18.5, 8;
     18, 9]);

% Define polytopes for the bottom fly zones
pursuer3_cvx1 = Polyhedron('V', ...
    [35.025, 1.25;
     35,-0.025;          % Central and Girard
     22.5, -0.0655;
     22.5, 1.25]);
pursuer3_cvx2 = Polyhedron('V', ...
    [25, 1.25;
     30, 1.25;
     30, 7.3;
     25.1, 7.3]);
pursuer3_cvx3 = Polyhedron('V', ...
    [25.5, 6.15;
     25.5, 9;
     20.5, 6.15;
     21.5, 9]);
 
% Plot
plot(area_of_interest, 'color', 'b', 'alpha',0.2);
hold on;
plot(pursuer1_cvx1, 'color', 'g', 'alpha',0.4);
plot(pursuer1_cvx2, 'color', 'g', 'alpha',0.4);
plot(pursuer2_cvx1, 'color', 'g', 'alpha',0.4);
plot(pursuer2_cvx2, 'color', 'g', 'alpha',0.4);
plot(pursuer2_cvx3, 'color', 'g', 'alpha',0.4);
plot(pursuer3_cvx1, 'color', 'g', 'alpha',0.4);
plot(pursuer3_cvx2, 'color', 'g', 'alpha',0.4);
plot(pursuer3_cvx3, 'color', 'g', 'alpha',0.4);
axis equal;
box on;
set(gca,'FontSize',20);
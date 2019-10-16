clf

% Define polytope for area of interest
area_of_interest = Polyhedron('V', ...
    [0, 0;              % University and Central
     5, -0.125;         % Terrace and Central
     35,-0.25;          % Central and Girard
     35.2, 2;           % Campus and Girard
     20.5, 5.5;         % Las Lomas and Campus
     1, 5]);            % Las Lomas and University

% Define polytopes for the left fly zones
polytope_flyzone_left_cvx1 = Polyhedron('V', ...
    [0, 0;              % University and Central
     1, 5;              % Las Lomas and University
     2.6, 5.04;         % Redondo and Las Lomas
     1.7, 0.5]);        % Bend of Redondo
polytope_flyzone_left_cvx2 = Polyhedron('V', ...
    [0, 0;              % University and Central
     1.7,  0.5;         % Bend of Redondo
     15,  0.25;         % Harvard into the University
     15,  -0.166;       % Harvard and Central
     5, -0.125]);       % Terrace and Central

% Define polytopes for the top fly zones
polytope_flyzone_top_cvx1 = Polyhedron('V', ...
    [8, 5.18;
     15, 5.358;
     15.1, 3.8;
     8, 3.6]);
polytope_flyzone_top_cvx2 = Polyhedron('V', ...
    [8, 3.6;
     7.25,3;
     13.65, 2.75;
     13.5, 3.75]);
polytope_flyzone_top_cvx3 = Polyhedron('V', ...
    [13.25, 3;
     14.5, 2.5;
     18.5, 2.5;
     18, 3]);

% Define polytopes for the bottom fly zones
polytope_flyzone_bottom_cvx1 = Polyhedron('V', ...
    [35.025, 0.25;
     35,-0.25;          % Central and Girard
     22.5, -0.2;
     22.5, 0.25]);
polytope_flyzone_bottom_cvx2 = Polyhedron('V', ...
    [25, 0.25;
     30, 0.25;
     30, 2.3;
     25.1, 2.3]);
polytope_flyzone_bottom_cvx3 = Polyhedron('V', ...
    [25.5, 2.15;
     25.5, 3;
     20.5, 2.15;
     21.5, 3]);
 
% Plot
plot(area_of_interest, 'color', 'b', 'alpha',0.2);
hold on;
plot(polytope_flyzone_left_cvx1, 'color', 'g', 'alpha',0.4);
plot(polytope_flyzone_left_cvx2, 'color', 'g', 'alpha',0.4);
plot(polytope_flyzone_top_cvx1, 'color', 'g', 'alpha',0.4);
plot(polytope_flyzone_top_cvx2, 'color', 'g', 'alpha',0.4);
plot(polytope_flyzone_top_cvx3, 'color', 'g', 'alpha',0.4);
plot(polytope_flyzone_bottom_cvx1, 'color', 'g', 'alpha',0.4);
plot(polytope_flyzone_bottom_cvx2, 'color', 'g', 'alpha',0.4);
plot(polytope_flyzone_bottom_cvx3, 'color', 'g', 'alpha',0.4);
axis square;
box on;
set(gca,'FontSize',20);

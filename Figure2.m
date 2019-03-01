clear
clc
close all
addpath('helperFuns/');

figure_number = 3; %Use 1-3

% Monte carlo parameters
no_of_grid_points = 1000;
no_of_MC_points = 1e5;
% Use this sigma for figures 1--3. Else generate new for figure 4
obstacle_sigma_temp = [11.6218    0.5928;
                       0.5928    3.7514];
fprintf('Generating Figure 2.%c\n',figure_number+64);
if figure_number == 1
    sampling_circle_radius = 50;
    obstacle_size = 25;
    obstacle_sigma_scaling = 1;
    obstacle_sigma = obstacle_sigma_scaling * obstacle_sigma_temp;
elseif figure_number == 2
    sampling_circle_radius = 25;
    obstacle_size = 10;
    obstacle_sigma_scaling = 1;
    obstacle_sigma = obstacle_sigma_scaling * obstacle_sigma_temp;
elseif figure_number == 3
    sampling_circle_radius = 35;
    obstacle_size = 10;
    obstacle_sigma_scaling = 2;
    obstacle_sigma = obstacle_sigma_scaling * obstacle_sigma_temp;
else
    sampling_circle_radius = 50;
    obstacle_size = 25;
    obstacle_sigma_scaling = 1;
    %% Obstacle sigma definition (randomly generated)
    obstacle_sigma_offset = [5,0;
                             0,1]; % Something greater than I_2
    obstacle_sigma_temp = rand(2,2);
    obstacle_sigma = 1 / 2 * ((obstacle_sigma_temp +...
                                    obstacle_sigma_temp') +...
                                  2* obstacle_sigma_offset);
end
alpha_value = 0.001;
no_of_points_on_the_circle = 10;
no_of_vectors_MinkSum = 50;
obstacle_mu = [2;
               2];

% Central symmetry being used in PrOccupySetConstraint function
obstacle_lower_bound = obstacle_size * [-1; 
                                        -1];
obstacle_upper_bound = -obstacle_lower_bound;

% Obstacle polytope and volume
obstacle_polytope = Polyhedron('lb', obstacle_lower_bound,...
                               'ub', obstacle_upper_bound);
obstacle_volume = volume(obstacle_polytope);

% Does the alpha pass the non-empty PrOccupySet test?
% \phi(x_{max}) > \alpha?
maximum_alpha = mvncdf(obstacle_lower_bound'+obstacle_mu',...
                       obstacle_upper_bound'+obstacle_mu',...
                       obstacle_mu',...
                       obstacle_sigma);
assert(alpha_value < maximum_alpha, 'Too high alpha value')

% fmincon settings
fminconsettings=optimset(@fmincon);
fminconsettings.Display='off';
%fminconsettings.Algorithm='sqp';
%fminconsettings.TolCon=1e-8;
%fminconsettings.TolFun=1e-8;

% Generating the Gaussian_level_set for visualization
disp('Generating the Gaussian_level_set for visualization')
[~,~,Gaussian_level_set] = ...
   getMinkSumBasedOverapproximationPrOccupySet(...
                                          obstacle_mu,...
                                          obstacle_sigma,...
                                          obstacle_polytope,...
                                          obstacle_volume,...
                                          sampling_circle_radius,...
                                          alpha_value,...
                                          no_of_vectors_MinkSum);

disp('Generating the MinkSumBased_overapproximation')
timer=tic;
%[MinkSumBased_overapproximation,...
 %Mink_boundary_point_of_interest] = ...
   %getMinkSumBasedOverapproximationPrOccupySet(...
                                          %obstacle_mu,...
                                          %obstacle_sigma,...
                                          %obstacle_polytope,...
                                          %obstacle_volume,...
                                          %sampling_circle_radius,...
                                          %alpha_value,...
                                          %no_of_points_on_the_circle);
[MinkSumSupportFunBased_overapproximation] = ...
getMinkSumSupportFunBasedOverapproximationPrOccupySetBox(...
                                       obstacle_mu,...
                                       obstacle_sigma,...
                                       obstacle_polytope,...
                                       alpha_value,...
                                       obstacle_volume,...
                                       20*no_of_points_on_the_circle);
elapsed_time_minksum = toc(timer);
disp('Generating the polytopic_overapproximation')
timer=tic;
[tight_polytope,...
 underapprox_polytope,...
 tight_boundary_point_of_interest] = getTightPolytopicOverapproximationPrOccupySet(...
                                                    obstacle_mu,...
                                                    obstacle_sigma,...
                                                    obstacle_lower_bound,...
                                                    obstacle_upper_bound,...
                                                    sampling_circle_radius,...
                                                    alpha_value,...
                                                    no_of_points_on_the_circle,...
                                                    fminconsettings);

elapsed_time_polytope = toc(timer);
fprintf('Is the overapproximation correct? ')
if tight_polytope.contains(underapprox_polytope)
    disp('Yes')
else
    disp('No')
end
% Should be zero since none of them should have make an acute angle with the projection line
%disp('Any acute angles?');
%disp(find(sum(tight_polytope.H*[tight_boundary_point_of_interest; -ones(1,no_of_points_on_the_circle)]>0)==1));
%% Should be all zeros since none of them should have a slack
%disp('Any slacked projections?');
%disp(find(diag(tight_polytope.H*[tight_boundary_point_of_interest; -ones(1,no_of_points_on_the_circle)])<0 == 1)');

timer=tic;
[xvec,yvec,frequency,other_values]=get_Gaussian_PrOccupySet_via_MonteCarlo(...
                                      no_of_MC_points,...
                                      no_of_grid_points,...
                                      obstacle_mu,...
                                      obstacle_sigma,...
                                      obstacle_size);
elapsed_time_MC = toc(timer);

%% Visual comparison
figure(1)
clf
plot(MinkSumSupportFunBased_overapproximation,'color','y','alpha',0.1);
hold on
plot(tight_polytope,'alpha',0.4,'color','g')
scatter(tight_boundary_point_of_interest(1,:),tight_boundary_point_of_interest(2,:),60, 'kd','filled');
plot(underapprox_polytope,'alpha',0.4,'color','b');
plot(Gaussian_level_set, 'alpha',0.8, 'color', 'k');
%plot(MinkSumBased_overapproximation, 'alpha',0.8, 'color', 'y');
%scatter(Mink_boundary_point_of_interest(1,:),Mink_boundary_point_of_interest(2,:), 60, 'rd', 'filled');
colormap('hot')
contour(xvec(1:10:end),yvec(1:10:end),frequency(1:10:end,1:10:end),[alpha_value alpha_value],'LineWidth',2);
leg = legend('$\mathrm{OvOccupySet}_{x}^+(\alpha)$',...
             '$\mathrm{OvOccupySet}_{x}(\alpha)$',...
             'Projection points on $\mathrm{OccupySet}_{x}(\alpha)$',...
             '$\mathrm{UnOccupySet}_{x}(\alpha)$',...
             '$\left\{ \overline{z}\in \mathcal{X}: \psi_{x}(\overline{x}) \geq \frac{\alpha}{ \mathrm{m}( \mathcal{G}_x(\overline{0}))}\right\}$',...
             'Monte-Carlo simulation'); %'Boundary points for MinkSum',...
%for point_index=1:no_of_points_on_the_circle
%    line_points = [faraway_points_backward(:,point_index),boundary_point_of_interest(:,point_index), faraway_points_forward(:,point_index)];
%    plot(line_points(1,:),line_points(2,:),'k-','LineWidth',1)
%end
set(leg,'Location','bestoutside','interpreter','latex')
grid on
box on
axis square
axis([-40 45 -35 40])
set(gca,'FontSize',20);
% xlabel('x')
% ylabel('y')

disp(elapsed_time_polytope)
disp(elapsed_time_minksum)
disp(elapsed_time_MC)

%% Ellipsoid is x^\top ellipsoid_shape^{-1} x <= 1
%probability_threshold = alpha_value/volume(obstacle_polytope);
%R_squared = - 2 *  log(sqrt(det(2*pi * obstacle_sigma)) * probability_threshold);
%Gaussian_level_set = ellipsoid(obstacle_mu, obstacle_sigma * R_squared);
%Gaussian_99_confidence = ellipsoid(obstacle_mu, obstacle_sigma * chi2inv(0.99,2));
%elloptions.color = 'g';
%elloptions.fill = 1;
%elloptions.width = 0.1;
%elloptions.shade = 0;
%plot(Gaussian_level_set, elloptions)

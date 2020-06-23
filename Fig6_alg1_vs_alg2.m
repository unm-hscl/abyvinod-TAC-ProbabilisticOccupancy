%% Use numbers {1,2,3} to obtain the different subfigures
% Use a fixed obstacle_sigma_temp (H in the text) for figures 1--3. Else,
% generate new if 4 is selected
% figure_number = 3; % Use 1-3

% Monte carlo parameters
no_of_grid_points = 1000;
no_of_MC_points = 1e3;
no_of_MC_points_high_acc = 1e5;
% Problem parameters
alpha_value = 0.001;
n_points_Alg1 = 10;
n_dir_vecs_Alg2 = 200;
obstacle_mu = [2; 2];
obstacle_sigma_temp = [11.62    0.59;
                        0.59    3.75];

fprintf('Generating Figure 6.%c\n',figure_number+64);
if figure_number == 1
    sampling_circle_radius = 50;
    obstacle_box_side_half_length = 25;
    obstacle_sigma_scaling = 1;
    obstacle_sigma = obstacle_sigma_scaling * obstacle_sigma_temp;
elseif figure_number == 2
    sampling_circle_radius = 25;
    obstacle_box_side_half_length = 10;
    obstacle_sigma_scaling = 1;
    obstacle_sigma = obstacle_sigma_scaling * obstacle_sigma_temp;
elseif figure_number == 3
    sampling_circle_radius = 35;
    obstacle_box_side_half_length = 10;
    obstacle_sigma_scaling = 2;
    obstacle_sigma = obstacle_sigma_scaling * obstacle_sigma_temp;
else
    % This block is never used! Just for testing and producing
    % obstacle_sigma_temp
    sampling_circle_radius = 50;
    obstacle_box_side_half_length = 25;
    obstacle_sigma_scaling = 1;
    %% Obstacle sigma definition (randomly generated)
    obstacle_sigma_offset = [5,0;
                             0,1]; % Something greater than I_2
    obstacle_sigma_temp = rand(2,2);
    obstacle_sigma = 1 / 2 * ((obstacle_sigma_temp +...
                                    obstacle_sigma_temp') +...
                                  2* obstacle_sigma_offset);
end

% Central symmetry being used in PrOccupySetConstraint function
obstacle_lower_bound = obstacle_box_side_half_length * [-1; -1];
obstacle_upper_bound = -obstacle_lower_bound;

% Obstacle polytope and volume
obstacle_polytope = Polyhedron('lb', obstacle_lower_bound,...
                               'ub', obstacle_upper_bound);
obstacle_polytope.minVRep();
obstacle_volume = volume(obstacle_polytope);

% Does the alpha pass the non-empty PrOccupySet test?
% \phi(x_{max}) > \alpha? --- x_max = obstacle_mu due to central symmetry
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

% Implement Algorithm 1
disp('Generating the polytopic_overapproximation')
timer=tic;
[tight_polytope, underapprox_polytope, tight_boundary_point_of_interest] = ...
    compAlg1Alg2_getTightPolytopicOverapproximationPrOccupySet(...
                                                    obstacle_mu,...
                                                    obstacle_sigma,...
                                                    obstacle_lower_bound,...
                                                    obstacle_upper_bound,...
                                                    sampling_circle_radius,...
                                                    alpha_value,...
                                                    n_points_Alg1,...
                                                    fminconsettings);
elapsed_time_polytope = toc(timer);
fprintf('Correct overapproximation and underapproximation returned by Alg. 1? ')
if tight_polytope.contains(underapprox_polytope)
    disp('Yes')
else
    disp('No')
    warning('Error!');
end
% Should be zero since none of them should have make an acute angle with the projection line
%disp('Any acute angles?');
%disp(find(sum(tight_polytope.H*[tight_boundary_point_of_interest; -ones(1,no_of_points_on_the_circle)]>0)==1));
%% Should be all zeros since none of them should have a slack
%disp('Any slacked projections?');
%disp(find(diag(tight_polytope.H*[tight_boundary_point_of_interest; -ones(1,no_of_points_on_the_circle)])<0 == 1)');

% Implement Algorithm 2
disp('Generating the MinkSumBased_overapproximation')
timer=tic;
[MinkSumSupportFunBased_overapproximation] = ...
    compAlg1Alg2_getMinkSumPrOccupySetBox(obstacle_mu,...
                                          obstacle_sigma,...
                                          obstacle_polytope,...
                                          alpha_value,...
                                          obstacle_volume,...
                                          n_dir_vecs_Alg2);
elapsed_time_minksum = toc(timer);


% Perform Monte-Carlo simulations
timer=tic;
[xvec,yvec,frequency,~] = ...
    compAlg1Alg2_getGaussianPrOccupySetMonteCarlo(no_of_MC_points,...
                                      no_of_grid_points,...
                                      obstacle_mu,...
                                      obstacle_sigma,...
                                      obstacle_box_side_half_length);
elapsed_time_MC = toc(timer);

timer=tic;
[~,~,frequency_high_acc,~] = ...
    compAlg1Alg2_getGaussianPrOccupySetMonteCarlo(no_of_MC_points_high_acc,...
                                      no_of_grid_points,...
                                      obstacle_mu,...
                                      obstacle_sigma,...
                                      obstacle_box_side_half_length);
elapsed_time_MC_high_acc = toc(timer);

%% Visual comparison
color_Mink = [52, 209, 191]/255;
color_tight = [239, 239, 239]/255;
color_under = [209, 52, 91]/255;
color_MC_low_acc = [47, 82, 224]/255;%/255;%[0.9290 0.6940 0.1250];
color_MC_high_acc = [0.4,0.4,0.4];%[14, 14, 14]/255;
figure(1)
clf
plot(MinkSumSupportFunBased_overapproximation,'color', color_Mink,'alpha', 1);
hold on
plot(tight_polytope,'alpha',1,'color', color_tight)
plot(underapprox_polytope,'alpha',1,'color', color_under);
plot(Polyhedron('lb',[-1,-1],'ub',[1,1]), 'alpha',1, 'color', 'k');   % Dummy polytope for Gaussian level set
%plot(MinkSumBased_overapproximation, 'alpha',0.8, 'color', 'y');
%scatter(Mink_boundary_point_of_interest(1,:),Mink_boundary_point_of_interest(2,:), 60, 'rd', 'filled');
contour(xvec(1:10:end),yvec(1:10:end),frequency(1:10:end,1:10:end),[alpha_value alpha_value], 'LineColor', color_MC_low_acc, 'LineWidth', 4);
contour(xvec(1:10:end),yvec(1:10:end),frequency_high_acc(1:10:end,1:10:end),[alpha_value alpha_value], 'LineColor', color_MC_high_acc, 'LineWidth', 4,'LineStyle','-.');
plot(Polyhedron('lb',[-1,-1],'ub',[1,1]), 'alpha',1, 'color', 'y');   % Dummy polytope for markers
leg = legend('$\mathrm{OccupySet}_{x}^{+\sharp}(\alpha)$',...
             '$\mathrm{OccupySet}_{x}^{\sharp}(\alpha)$',...
             '$\mathrm{OccupySet}_{x}^{\flat}(\alpha)$',...
             '$\left\{ \overline{x}\in \mathcal{X}: \psi_{x}(\overline{x}) \geq \frac{\alpha}{ \mathrm{m}( \mathcal{G}_x)}\right\}$',...
             'Sampling ($10^3$ samples)', ...
             'Sampling ($10^5$ samples)', ...
             'Projection points'); %'Boundary points for MinkSum',...
%for point_index=1:no_of_points_on_the_circle
%    line_points = [faraway_points_backward(:,point_index),boundary_point_of_interest(:,point_index), faraway_points_forward(:,point_index)];
%    plot(line_points(1,:),line_points(2,:),'k-','LineWidth',1)
%end
set(leg,'Location','bestoutside','interpreter','latex','AutoUpdate','off')
scatter(tight_boundary_point_of_interest(1,:),tight_boundary_point_of_interest(2,:), 150, 'ys','filled','MarkerEdgeColor','k');
% Draw the level set
level_set_threshold = alpha_value/obstacle_volume;
compAlg1Alg2_plotGaussianLevelSet(obstacle_mu, obstacle_sigma, ...
    level_set_threshold)
grid on
box on
axis square
axis([-40 45 -35 40])
% Zoom the figure size till we have the legend occupy 3/4 of y-axis
set(gca,'FontSize', 45);
xlabel('x')
ylabel('y')
disp('Computation time in seconds');
disp('Truth |  MC    | Alg. 1| Alg. 2');
fprintf('%2.3f | %2.3f | %2.3f | %2.3f\n', elapsed_time_MC_high_acc, elapsed_time_MC, elapsed_time_polytope, elapsed_time_minksum);
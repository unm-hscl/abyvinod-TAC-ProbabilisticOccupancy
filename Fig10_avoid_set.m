clear
close all;
clc;
addpath('helperFuns/Fig10/');

figure_number = 4; %3-4
fontSize = 30;
sampling_time = 0.05;
probability_threshold = 1 - 0.99;
no_of_vectors_MinkSum = 50;

%% Input velocity --- Gaussian                      
mu_velocity = 5;
var_velocity = 1;
v_random_vector = RandomVector.gaussian(mu_velocity, var_velocity);

%% Obstacle geometry
obstacle_radius = 0.5;
obstacle_polytope = Polyhedron('V', obstacle_radius * [cos(linspace(0,2*pi,50))',sin(linspace(0,2*pi,50))']);
obstacle_volume = pi * obstacle_radius^2;

%% Obstacle definition
obstacle_init_location = [10;
                          10];
omega_init = 0;
theta_init = pi/4;

%% Turning rate stochasticity
omega_min = -20;
% Q is symmetric about origin
omega_max = - omega_min;
% Create omega values possible (mode set Q) and its |.|
omega_values = linspace(omega_min, omega_max, 5);
omega_step = omega_values(2) - omega_values(1);
no_of_omega_values = length(omega_values);
% Compute the index of omega_init
indx_of_omega_init = round((omega_init - omega_min)/omega_step) + 1;
% Determine the pdf governing the turning rate
if figure_number == 4
    disp('Figure 3.b data used!');
    % Skewed towards anti-clockwise turning
    omega_pdf = [0.5 0.45 0.05 0.0 0.0];
else
    disp('Figure 3.a data used!');
    % Equally spread out
    omega_pdf = 1/no_of_omega_values*ones(1,no_of_omega_values);    
end


%% Mode stays on for \tausw time steps
tausw = 5;
last_time_step_including_init = 15;       % Last time step
% For time_horizon = 3, time_steps_taken = 4 (includes t=0)
%no_of_switches = floor((last_time_step - 1)/tausw);
%% Generate all possible combination switches --- TODO: Shift to iterator
% All combinations from the first switch onwards
all_comb_indx_row_wise_temp = allcomb(1:no_of_omega_values,1:no_of_omega_values);
no_of_sequences = size(all_comb_indx_row_wise_temp,1);
% Prepend the initial mode state
all_comb_indx_row_wise = [repmat(indx_of_omega_init, no_of_sequences, 1),...
                          all_comb_indx_row_wise_temp];

% Translated to turning rate sequences
% Use transpose since omega_values is horizontal
all_comb_omega_row_wise =...
                [omega_values(all_comb_indx_row_wise(:,1))', ...
                 zeros(no_of_sequences, tausw -1), ...
                 omega_values(all_comb_indx_row_wise(:,2))', ...
                 zeros(no_of_sequences, tausw -1), ...
                 omega_values(all_comb_indx_row_wise(:,3))', ...
                 zeros(no_of_sequences, tausw -1)];  
% Translated to heading sequences
diff_heading = cumsum(all_comb_omega_row_wise * sampling_time,2);
all_comb_heading_row_wise = repmat(theta_init, no_of_sequences, last_time_step_including_init) + ...
    diff_heading;

%% Get the probability for the particular mode sequence
% Probability that omegas at the switching instants are exactly what we
% specified
probabilities_of_sequences = prod(omega_pdf(all_comb_indx_row_wise(:,2:3)),2);
nnz_prob_sequences = find(probabilities_of_sequences > 0);
no_nnz_prob_sequences = length(nnz_prob_sequences);

% Additional definitions
probability_occupied_set_overapprox = cell(no_of_sequences,1);
for indx = 1:no_of_sequences
    probability_occupied_set_overapprox{indx} = Polyhedron();
end
elapsed_time = 0;

%% Iterate over all non-zero probability sequences
disp('Indx | Mode_sequence |  Mode_Prob  | Status');    
for indx_indx_sequence = 1:no_nnz_prob_sequences
    % Translate the indx to indx_sequence
    indx_sequence = nnz_prob_sequences(indx_indx_sequence);
    % Create obstacle dynamics associated with this theta sequence
    theta_seq = all_comb_heading_row_wise(indx_sequence,:);
    sys = LtvSystem('StateMatrix', @(t) eye(2), ...
                    'DisturbanceMatrix', @(t) sampling_time * ...
                        [cos(theta_seq(t+1)); sin(theta_seq(t+1))], ...
                    'Disturbance', v_random_vector);
    fprintf(' %2d  | %13s |   %1.4f    | \n', indx_sequence, ...
        num2str(all_comb_indx_row_wise(indx_sequence,:)), ...
        probabilities_of_sequences(indx_sequence));
    % Compute
    timerVal = tic;
    % Compute alpha_sq
    alpha_q_tau = probability_threshold /...
        (probabilities_of_sequences(indx_sequence) * no_nnz_prob_sequences);
    if alpha_q_tau > 1
        % The PrOccupySet for the corresponding DTPV is empty since
        % probabilistic occupancy function = 1 is also not enough            
        probability_occupied_set_overapprox{indx_sequence} = Polyhedron();
        fprintf('\b Skipped as alpha_S >1!\n')
    else
        % Get the mean and covariance of the state at time tau
        state_at_tau = SReachFwd('state-stoch', sys, obstacle_init_location, ...
            last_time_step_including_init);
        obstacle_mu = state_at_tau.mean();
        obstacle_sigma = state_at_tau.cov();
        % Compute the maxima of the occupancy function (Overapproximated by a
        % circumscribing box)
        if det(obstacle_sigma)>1e-8
            mode_of_occupancy_function = mvncdf(-obstacle_radius + obstacle_mu',...
                                                 obstacle_radius + obstacle_mu',...
                                                obstacle_mu',...
                                                obstacle_sigma);
        else
            % Specialized for 2D
            [V,D] = eig(obstacle_sigma);
            direction_of_sigma = V(:,2);
            sigma_value = sqrt(D(2,2));
            mode_of_occupancy_function = normcdf(obstacle_radius,0,sigma_value)-normcdf(-obstacle_radius,0,sigma_value);
        end
        if mode_of_occupancy_function < alpha_q_tau
            % Avoid set is empty! Don't plot the trajectory
            probability_occupied_set_overapprox{indx_sequence} = Polyhedron();
            fprintf('\b Skip it! Mode: %1.3f < alpha_S: %1.3f\n', mode_of_occupancy_function, alpha_q_tau);
        else
            fprintf('\b Compute! Mode: %1.3f > alpha_S: %1.3f\n', mode_of_occupancy_function, alpha_q_tau)
            % Compute the probability_occupied_set
            MinkSumSupportFunBased_overapproximation =...
                getMinkSumSupportFunBasedOverapproximationPrOccupySet(...
                                                obstacle_mu,...
                                                obstacle_sigma,...
                                                obstacle_radius,...
                                                alpha_q_tau,...
                                                obstacle_volume,...
                                                no_of_vectors_MinkSum);
            % Add the non-empty probability_occupied_set to the union
            probability_occupied_set_overapprox{indx_sequence} =...
                            MinkSumSupportFunBased_overapproximation;
        end  
    end
    elapsed_time = elapsed_time + toc(timerVal);
end

figure(1)
clf
axis square
%set(gca,'XTick',5:2:15)
%set(gca,'YTick',5:2:15)
box on
grid on
xlabel('x');
ylabel('y');
set(gca,'FontSize', fontSize);
color_scatter = ['bx';'bo';'bd';'bs';'b+';'rx';'ro';'rd';'rs';'r+';'kx';'ko';'kd';'ks';'k+';'mx';'mo';'md';'ms';'m+';'cx';'co';'cd';'cs';'c+';'yx';'yo';'yd';'ys';'y+';'gx';'go';'gd';'gs';'g+';];
hold on
for indx_indx_sequence = 1:no_nnz_prob_sequences
    % Translate the indx to indx_sequence
    indx_sequence = nnz_prob_sequences(indx_indx_sequence);
    % Create obstacle dynamics associated with this theta sequence
    theta_seq = all_comb_heading_row_wise(indx_sequence,:);
    sys = LtvSystem('StateMatrix', @(t) eye(2), ...
                    'DisturbanceMatrix', @(t) sampling_time * ...
                        [cos(theta_seq(t+1)); sin(theta_seq(t+1))], ...
                    'Disturbance', v_random_vector);
    % Get the polytope
    proccupyPolytope = probability_occupied_set_overapprox{indx_sequence};
    % Plot only if it is non-empty
    if ~isEmptySet(proccupyPolytope)
        plot(proccupyPolytope,'color',color_scatter(indx_sequence,1),'alpha',0.1);
        scatter(obstacle_init_location(1), obstacle_init_location(2),color_scatter(indx_sequence,1:2));                
        for time_in_evolution=1:last_time_step_including_init
            % Get the mean and covariance of the state at time tau
            state_at_tau = SReachFwd('state-stoch', sys, obstacle_init_location, ...
                time_in_evolution);
            mu_point = state_at_tau.mean();
            scatter(mu_point(1),mu_point(2),color_scatter(indx_sequence,1:2));                
        end
    end
end
checkAvoidSetUsingMonteCarlo
%% Report times
fprintf('Elapsed time\n    MC | Alg. 2\n %1.3f | %1.3f\n', elapsed_time_MC, ...
    elapsed_time)
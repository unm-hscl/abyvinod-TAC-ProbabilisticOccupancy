clear
clc

timer=tic;
sampling_time = 0.05;
probability_threshold = 1 - 0.99;
no_of_points_on_the_circle = 10;

% fmincon settings
fminconsettings=optimset(@fmincon);
fminconsettings.Display='off';
%fminconsettings.Algorithm='sqp';
%fminconsettings.TolCon=1e-8;
%fminconsettings.TolFun=1e-8;


% Additional definitions
probabilities_of_sequences=[];      
probability_occupied_set_overapprox_union = PolyUnion();
probability_occupied_set_underapprox_union = PolyUnion();
total_time_without_plotting = 0;

%% Obstacle definition
obstacle_init_location = [10;
                          10];
mu_velocity = 5;
variance_velocity = 1;
omega_min = -10;
omega_step = 5;
omega_init = 0;
obstacle_size = 1;
sampling_circle_radius = 100;
obstacle_lower_bound = obstacle_size * [-1; 
                                        -1];
theta_init = pi/4;
% Mode stays on for time steps 
% how many steps does sequence mask stay on (inclusive of t=0)
sequence_mask = [3 3 3];            
% For time_horizon = 3, last_count_of_time_step = 4 (includes t=0)
last_count_of_time_step = sum(sequence_mask); 
% Q is symmetric about origin
omega_max = - omega_min;
% Central symmetry being used in PrOccupySetConstraint function
obstacle_upper_bound = -obstacle_lower_bound;
% Obstacle polytope and volume
obstacle_polytope = Polyhedron('lb', obstacle_lower_bound,...
                               'ub', obstacle_upper_bound);
obstacle_volume = volume(obstacle_polytope);
% Create omega values possible (mode set Q) and its |.|
omega_values = omega_min: omega_step : -omega_min;
no_of_omega_values = length(omega_values);
% Compute the index of omega_init
indx_of_omega_init = round((omega_init - omega_min)/omega_step) + 1;
% Markov transition matrix each point has equal likelihood
 switching_probability_matrix = 1/no_of_omega_values*ones(no_of_omega_values, no_of_omega_values);
% % Markov transition matrix favouring turning -\omega the most
%switching_probability_matrix = [0.5 0.2 0.1 0.1 0.1;
                                %0.5 0.2 0.1 0.1 0.1;
                                %0.5 0.2 0.1 0.1 0.1;
                                %0.5 0.2 0.1 0.1 0.1;
                                %0.5 0.2 0.1 0.1 0.1];
%switching_probability_matrix = [0.5 0.45 0.04 0.005 0.005;
                                %0.5 0.45 0.04 0.005 0.005;
                                %0.5 0.45 0.04 0.005 0.005;
                                %0.5 0.45 0.04 0.005 0.005;
                                %0.5 0.45 0.04 0.005 0.005];
% % Random Markov transition matrix
% switching_probability_matrix = randfixedsum(no_of_omega_values, no_of_omega_values, 1, 0, 1)';
% Check if the stochastic matrix has the correct structure (rows sum up to 1)
rowsum_of_matrix = sum(switching_probability_matrix,2);
assert(sum(abs(rowsum_of_matrix - 1))<1e-8,...
                'Socbox:InvalidArgs',...
                'Invalid stochastic matrix: Rows do not sum to one.');
assert(size(switching_probability_matrix,1) == no_of_omega_values,...
                'Socbox:InvalidArgs',...
                'Invalid stochastic matrix: Size is incorrect.')

%% Generate all possible combination switches --- TODO: Shift to iterator
% All combinations from the first switch onwards
all_comb_indx_row_wise_temp = allcomb(1:no_of_omega_values,1:no_of_omega_values);
no_of_sequences = size(all_comb_indx_row_wise_temp,1);
% Prepend the initial mode state
all_comb_indx_row_wise = [repmat(indx_of_omega_init, no_of_sequences, 1),...
                          all_comb_indx_row_wise_temp];
assert(sequence_mask(1)>1,...
       'Socbox:InvalidArgs',...
       'Initial omega stays on at least once (t=0)');
% Translated to mode
all_comb_mode_row_wise =...
                [repmat(omega_init, no_of_sequences, sequence_mask(1)-1),...
                 repmat(omega_values(all_comb_indx_row_wise(:,2))',1,sequence_mask(2)),...
                 repmat(omega_values(all_comb_indx_row_wise(:,3))',1,sequence_mask(3))];

%% Iterate over all sequences => Compute their probability and probability_occupied_set (if enough weight)
disp('Indx | Mode_sequence | Mode_Prob | Status');    
for indx_sequence = 1: no_of_sequences
    timer = tic;
    % Get the probability for the particular mode sequence
    probabilities_of_sequences(indx_sequence) =...
                            get_probability_for_mode_sequence(...
                                    all_comb_indx_row_wise(indx_sequence,:),...
                                    switching_probability_matrix,...
                                    indx_of_omega_init);
    fprintf(' %2d  | %13s |   %1.2f    | \n',indx_sequence, num2str(all_comb_indx_row_wise(indx_sequence,:)),probabilities_of_sequences(indx_sequence));
    probability_threshold_for_sequence = probability_threshold /...
                                (probabilities_of_sequences(indx_sequence) * ...
                                    no_of_sequences);
    % Get the mode sequence
    mode_sequence_omega = all_comb_mode_row_wise(indx_sequence,:);
    if probability_threshold_for_sequence > 1
        % The PrOccupySet for the corresponding DTPV is empty since
        % probabilistic occupancy function = 1 is also not enough            
        total_time_without_plotting = total_time_without_plotting + toc(timer);
        fprintf('\b Skipped as alpha_S >1!\n')
    else
        % Compute the ctrb and state_transition_matrix for the underlying DTPV
        [ctrb_matrix, state_transition_matrix] =...
        get_ctrb_and_state_transition_matrices_unicycle(mode_sequence_omega,...
                                                        theta_init,...
                                                        sampling_time);
        % Compute mu and sigma for the obstacle FSRPD at the time instant
       [obstacle_mu, obstacle_sigma] = get_FSRPD_mean_and_covariance_matrix(...
                                            ctrb_matrix,...
                                            state_transition_matrix,...
                                            mu_velocity,...
                                            variance_velocity,...
                                            last_count_of_time_step,...
                                            obstacle_init_location);
        % Compute the maxima of the occupancy function
        mode_of_occupancy_function = mvncdf(obstacle_lower_bound'+obstacle_mu',...
                                            obstacle_upper_bound'+obstacle_mu',...
                                            obstacle_mu',...
                                            obstacle_sigma);
        if mode_of_occupancy_function < probability_threshold_for_sequence
            % Avoid set is empty! Don't plot the trajectory
            fprintf('\b Skipped Mode: %1.3f < alpha_S: %1.3f !\n', mode_of_occupancy_function, probability_threshold_for_sequence);
        else
            fprintf('\b Computed Mode: %1.3f > alpha_S: %1.3f !\n', mode_of_occupancy_function, probability_threshold_for_sequence)
            if indx_sequence ~= 13
                % Compute the probability_occupied_set
                [MinkSumBased_overapproximation,...
                 Mink_boundary_point_of_interest] = ...
                   getMinkSumBasedOverapproximationPrOccupySet(...
                                                          obstacle_mu,...
                                                          obstacle_sigma,...
                                                          obstacle_polytope,...
                                                          obstacle_volume,...
                                                          sampling_circle_radius,...
                                                          probability_threshold_for_sequence,...
                                                          no_of_points_on_the_circle);
                % Add the non-empty probability_occupied_set to the union
                probability_occupied_set_overapprox_union.add(MinkSumBased_overapproximation);

%                % Compute the probability_occupied_set
%                [tight_polytope,...
%                 underapprox_polytope,...
%                 tight_boundary_point_of_interest] = ...
%                        getTightPolytopicOverapproximationPrOccupySet(...
%                                                obstacle_mu,...
%                                                obstacle_sigma,...
%                                                obstacle_lower_bound,...
%                                                obstacle_upper_bound,...
%                                                sampling_circle_radius,...
%                                                probability_threshold_for_sequence,...
%                                                no_of_points_on_the_circle,...
%                                                fminconsettings);
%                % Add the non-empty probability_occupied_set to the union
%                probability_occupied_set_overapprox_union.add(tight_polytope);
%                probability_occupied_set_underapprox_union.add(underapprox_polytope);
            end
        end
    end
    total_time_without_plotting = total_time_without_plotting + toc(timer);
end
%            % Compute mu_trajectory from t=0 to time_horizon
%            for time_in_evolution=1:length(mode_sequence_omega)+1
%                % Plot the mean trajectory
%                scatter(mu_trajectory{time_in_evolution}(1),mu_trajectory{time_in_evolution}(2),color_scatter(indx_sequence,1:2));                
%            end
%            % Plot the decision epochs
%            for time_in_evolution = 1:time_steps_between_switches:length(mode_sequence_omega)+1
%                scatter(mu_trajectory{time_in_evolution}(1),mu_trajectory{time_in_evolution}(2),130,color_scatter(indx_sequence,1:2));                
%            end

figure(2)
plot(probability_occupied_set_underapprox_union,'alpha',0);
axis square
axis([5 15 5 15])
set(gca,'XTick',5:2:15)
set(gca,'YTick',5:2:15)


figure(1)
plot(probability_occupied_set_overapprox_union,'alpha',0);
axis square
axis([5 15 5 15])
set(gca,'XTick',5:2:15)
set(gca,'YTick',5:2:15)
box on

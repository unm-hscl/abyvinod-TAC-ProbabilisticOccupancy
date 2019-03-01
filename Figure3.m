clear
close all;
clc;
addpath('helperFuns/');

figure_number = 4; %3-4

sampling_time = 0.05;
probability_threshold = 1 - 0.99;
no_of_vectors_MinkSum = 50;

% fmincon settings
fminconsettings=optimset(@fmincon);
fminconsettings.Display='off';
%fminconsettings.Algorithm='sqp';
%fminconsettings.TolCon=1e-8;
%fminconsettings.TolFun=1e-8;


% Additional definitions
probabilities_of_sequences=[];      
probability_occupied_set_overapprox = {};
probability_occupied_set_underapprox = {};
elapsed_time = 0;

%% Obstacle definition
obstacle_init_location = [10;
                          10];
mu_velocity = 5;
variance_velocity = 1;
omega_min = -5;
omega_step = 2.5;
omega_init = 0;
obstacle_radius = 0.5;
obstacle_polytope = Polyhedron('V', obstacle_radius * [cos(linspace(0,2*pi,50))',sin(linspace(0,2*pi,50))']);
obstacle_volume = pi * obstacle_radius^2;
theta_init = pi/4;
% Mode stays on for time steps 
% how many steps does sequence mask stay on (inclusive of t=0)
sequence_mask = [5 5 5];            
% For time_horizon = 3, time_steps_taken = 4 (includes t=0)
time_steps_taken = sum(sequence_mask); 
% Q is symmetric about origin
omega_max = - omega_min;
% Create omega values possible (mode set Q) and its |.|
omega_values = omega_min: omega_step : -omega_min;
no_of_omega_values = length(omega_values);
% Compute the index of omega_init
indx_of_omega_init = round((omega_init - omega_min)/omega_step) + 1;
if figure_number == 4
    disp('Figure 3.b data used!');
    % % Markov transition matrix favouring turning -\omega the most
    %switching_probability_matrix = [0.35 0.0 0.3 0.0 0.35;
                                    %0.35 0.0 0.3 0.0 0.35;
                                    %0.35 0.0 0.3 0.0 0.35;
                                    %0.35 0.0 0.3 0.0 0.35;
                                    %0.35 0.0 0.3 0.0 0.35];
    switching_probability_matrix = [0.5 0.47 0.03 0.0 0.0;
                                    0.5 0.47 0.03 0.0 0.0;
                                    0.5 0.47 0.03 0.0 0.0;
                                    0.5 0.47 0.03 0.0 0.0;
                                    0.5 0.47 0.03 0.0 0.0];
else
    disp('Figure 3.a data used!');
    % Markov transition matrix each point has equal likelihood
    switching_probability_matrix = 1/no_of_omega_values*ones(no_of_omega_values, no_of_omega_values);
end
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
    timerVal = tic;
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
        probability_occupied_set_overapprox{indx_sequence} = Polyhedron();
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
                                            time_steps_taken,...
                                            obstacle_init_location);
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
        if mode_of_occupancy_function < probability_threshold_for_sequence
            % Avoid set is empty! Don't plot the trajectory
            probability_occupied_set_overapprox{indx_sequence} = Polyhedron();
            fprintf('\b Skip it! Mode: %1.3f < alpha_S: %1.3f\n', mode_of_occupancy_function, probability_threshold_for_sequence);
        else
            fprintf('\b Compute! Mode: %1.3f > alpha_S: %1.3f\n', mode_of_occupancy_function, probability_threshold_for_sequence)
            % Compute the probability_occupied_set
            MinkSumSupportFunBased_overapproximation =...
                getMinkSumSupportFunBasedOverapproximationPrOccupySet(...
                                                obstacle_mu,...
                                                obstacle_sigma,...
                                                obstacle_radius,...
                                                probability_threshold_for_sequence,...
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
set(gca,'FontSize',20);
color_scatter = ['bx';'bo';'bd';'bs';'b+';'rx';'ro';'rd';'rs';'r+';'kx';'ko';'kd';'ks';'k+';'mx';'mo';'md';'ms';'m+';'cx';'co';'cd';'cs';'c+';'yx';'yo';'yd';'ys';'y+';'gx';'go';'gd';'gs';'g+';];
hold on
for indx_sequence = 1: no_of_sequences %no_of_sequences %
    %if mod(indx_sequence,3)==1
    proccupyPolytope = probability_occupied_set_overapprox{indx_sequence};
    if ~isempty(proccupyPolytope) && ~isEmptySet(proccupyPolytope)
        plot(proccupyPolytope,'color',color_scatter(indx_sequence,1),'alpha',0.1);
        mode_sequence_omega = all_comb_mode_row_wise(indx_sequence,:);
        [ctrb_matrix, state_transition_matrix] =...
        get_ctrb_and_state_transition_matrices_unicycle(mode_sequence_omega,...
                                                        theta_init,...
                                                        sampling_time);
        scatter(obstacle_init_location(1), obstacle_init_location(2),color_scatter(indx_sequence,1:2));                
        for time_in_evolution=1:length(mode_sequence_omega)+1
            % state_transition_matrix is eye(2) anyways
            [mu_point, ~] = get_FSRPD_mean_and_covariance_matrix(...
                                                ctrb_matrix(:,1:time_in_evolution),...
                                                eye(2),...
                                                mu_velocity,...
                                                variance_velocity,...
                                                time_in_evolution,...
                                                obstacle_init_location);
            scatter(mu_point(1),mu_point(2),color_scatter(indx_sequence,1:2));                
        end
    end
end
checkAvoidSetUsingMonteCarlo
disp(elapsed_time)
disp(elapsed_time_MC)

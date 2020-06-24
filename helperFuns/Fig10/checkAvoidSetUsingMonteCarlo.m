% Reusing: 
% Variables:
% mode_sequences_omega, 
% mu_velocity,
% variance_velocity, 
% all_comb_mode_row_wise,
% probabilities_of_sequences,
% sampling_time,
% theta_init
% Functions:
% get_ctrb_and_state_transition_matrices_unicycle
timerVal=tic;
grid_x_space = 0.05;
grid_y_space = 0.05;
grid_min_x = 8;
grid_max_x = 15;
grid_min_y = 8;
grid_max_y = 15;
xvec = grid_min_x:grid_x_space:grid_max_x;
yvec = grid_min_y:grid_y_space:grid_max_y;
no_of_grid_points_x = length(xvec);
no_of_grid_points_y = length(yvec);

no_of_particles_as_per_probability_seq = round(probabilities_of_sequences*no_of_MC_particles);
no_of_MC_particles = sum(no_of_particles_as_per_probability_seq);
%% Create the MC particles
MC_counter = zeros(no_of_grid_points_y,no_of_grid_points_x);

%% Propagate the dynamics
mode_sequence_indx = 0;
fprintf('Simulating the trajectories\nStops at %d (no_of_particles allocated)\n', no_of_sequences);
for indx_indx_sequence = 1:no_nnz_prob_sequences
    % Translate the indx to indx_sequence
    indx_sequence = nnz_prob_sequences(indx_indx_sequence);
    % Create obstacle dynamics associated with this theta sequence
    theta_seq = all_comb_heading_row_wise(indx_sequence,:);
    sys = LtvSystem('StateMatrix', @(t) eye(2), ...
                    'DisturbanceMatrix', @(t) sampling_time * ...
                        [cos(theta_seq(t+1)); sin(theta_seq(t+1))], ...
                    'Disturbance', v_random_vector);
    no_of_particles = no_of_particles_as_per_probability_seq(indx_sequence);
    [Z, ~, G] = sys.getConcatMats(last_time_step_including_init);
    state_transition_matrix = Z(end-1:end,:);
    ctrb_matrix = G(end-1:end,:);
    fprintf(' %2d ( %6d)\n', indx_sequence, no_of_particles);
    velocity_random_variable_realizations = sqrt(var_velocity)*randn(last_time_step_including_init, no_of_particles) + mu_velocity;
    particle_terminal_states =...
            state_transition_matrix * repmat(obstacle_init_location,1,no_of_particles)...
            + ctrb_matrix * velocity_random_variable_realizations;
    for x_indx = 1:length(xvec)
        for y_indx = 1:length(yvec)
            delta_location = particle_terminal_states - repmat([xvec(x_indx);yvec(y_indx)], 1, no_of_particles);
            % Evaluate probability of x_\tau\in \mathcal{O{ point_of_interest})
            MC_counter(y_indx,x_indx) = MC_counter(y_indx,x_indx) + sum(sum(delta_location.^2,1)<= obstacle_radius^2);
        end
    end
    %scatter(particle_terminal_states(1,:),particle_terminal_states(2,:),30,'ks')
end
frequency = MC_counter/no_of_MC_particles;
elapsed_time_MC = toc(timerVal);
fprintf('\n');
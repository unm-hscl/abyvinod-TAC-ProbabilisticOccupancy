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
no_of_MC_particles = 1e5;
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
for indx_sequence=1:no_of_sequences
    mode_sequence_omega = all_comb_mode_row_wise(indx_sequence,:);
    [ctrb_matrix, state_transition_matrix] =...
    get_ctrb_and_state_transition_matrices_unicycle(mode_sequence_omega,...
                                                    theta_init,...
                                                    sampling_time);
    no_of_particles = no_of_particles_as_per_probability_seq(indx_sequence);
    fprintf(' %2d ( %6d)\n', indx_sequence, no_of_particles);
    velocity_random_variable_realizations = sqrt(variance_velocity)*randn(time_steps_taken,no_of_particles) + mu_velocity;
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

colormap([1 0 0]);
contour(xvec,yvec,frequency,[probability_threshold probability_threshold],'linewidth',5)
axis([grid_min_x grid_max_x grid_min_y grid_max_y])
figure(2);
clf
surf(xvec,yvec,frequency);
axis square
axis([grid_min_x grid_max_x grid_min_y grid_max_y])
%set(gca,'XTick',5:2:15)
%set(gca,'YTick',5:2:15)
box on
grid on
set(gca,'FontSize',20);

grid_x_space = 0.05;
grid_y_space = 0.05;
grid_min_x = 9;
grid_max_x = 14.5;
grid_min_y = 9;
grid_max_y = 14.5;
xvec = grid_min_x:grid_x_space:grid_max_x;
yvec = grid_min_y:grid_y_space:grid_max_y;
no_of_grid_points_x = length(xvec);
no_of_grid_points_y = length(yvec);

mu_collection = {};
sigma_collection = {};

for indx_sequence=1:no_of_sequences
    mode_sequence_omega = all_comb_mode_row_wise(indx_sequence,:);
    [ctrb_matrix, state_transition_matrix] =...
    get_ctrb_and_state_transition_matrices_unicycle(mode_sequence_omega,...
                                                    theta_init,...
                                                    sampling_time);
    % Compute mu and sigma for the obstacle FSRPD at the time instant
    [mu_collection{indx_sequence}, sigma_collection{indx_sequence}] =...
                    get_FSRPD_mean_and_covariance_matrix(...
                                        ctrb_matrix,...
                                        state_transition_matrix,...
                                        mu_velocity,...
                                        variance_velocity,...
                                        time_steps_taken,...
                                        obstacle_init_location);
     
end
occupancy_function_values = zeros(no_of_grid_points_y,no_of_grid_points_x);
for x_indx=1:no_of_grid_points_x
    for y=vec
        point_of_interest = [x;y];
        
    end
end

% Target dynamics
srtinit

target_init_state = [15;0;-0.5;0];
relv_states = [1,3];
n_monte_carlo = 1e3;
skip_mc = 1;                        %Skip MC particles
dist_peak = 0;
dist_delta = 2;
dist_min = dist_peak - dist_delta;
dist_max = dist_peak + dist_delta;
plot_t_skip = 3;

% time_horizon = 10;
% optimal_input_vec = 10*repmat([5;1], time_horizon, 1);
target_affine_vec = [repmat([-7;0.5], 10, 1);
                     repmat([15;15], 5, 1);
                     repmat([15;-0.5], 10, 1)];
time_horizon = length(target_affine_vec)/2;

% Target system definition
sys_DI_1D = getChainOfIntegLtiSystem(2, 0.1, Polyhedron('lb',-1,'ub',1));
target_sys_state_mat = blkdiag(sys_DI_1D.state_mat, sys_DI_1D.state_mat);
target_sys_input_mat = blkdiag(sys_DI_1D.input_mat, sys_DI_1D.input_mat);
target_sys_dist_mat = blkdiag(sys_DI_1D.input_mat, sys_DI_1D.input_mat);
pd = makedist('Triangular','a', dist_min,'b', dist_peak,'c', dist_max);
target_sys_dist = RandomVector('UserDefined', @(N) random(pd,2,N));
target_sys = LtiSystem('StateMatrix', target_sys_state_mat, ...
                       'InputMatrix', target_sys_input_mat, ...
                       'DisturbanceMatrix', target_sys_dist_mat, ...
                       'Disturbance', target_sys_dist, ...
                       'InputSpace', ...
                        Polyhedron('lb', -[inf,inf], 'ub', [inf,inf]));

% Generate Monte Carlo simulation                    
target_concat_state_realization = generateMonteCarloSims(n_monte_carlo, ...
    target_sys, target_init_state, time_horizon, target_affine_vec);

% Target system support computation
Figure4_target_support
for t_indx = 2:plot_t_skip:time_horizon+1
    % Time goes from 0 to time_horizon for both
    %   target_support_position and target_concat_state_realization
    fprintf('Plotting time: %d\n', t_indx-1);
    plot(target_support_position(t_indx), 'alpha', 0.2, 'color', 'y');
    relv_indx = 4*(t_indx-1) + relv_states;
    scatter(target_concat_state_realization(relv_indx(1),1:skip_mc:end), ...
            target_concat_state_realization(relv_indx(2),1:skip_mc:end), ...
            'ro', 'filled');
    drawnow
end

% Target system distribution computation
% [set_of_polytopes, hpoly] = polytopesFromMonteCarloSims( ...
%      concat_state_realization, 4, relv_states, ...
%      {'color', 'red','alpha',0.2}, plot_t_skip);

% Plot mean trajectory
mean_trajectory = mean(target_concat_state_realization,2);
plot(mean_trajectory(relv_states(1):4:end), mean_trajectory(relv_states(2):4:end), 'r--', 'linewidth', 2);
                            

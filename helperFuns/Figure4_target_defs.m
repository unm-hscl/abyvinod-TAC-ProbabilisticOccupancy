% Target dynamics
srtinit

target_init_state = [15;0;-0.5;0];
n_monte_carlo = 1e3;
skip_monte_carlo = 1;
dist_peak = 0;
dist_delta = 2;
dist_min = dist_peak - dist_delta;
dist_max = dist_peak + dist_delta;
plot_t_skip = 3;

% time_horizon = 10;
% optimal_input_vec = 10*repmat([5;1], time_horizon, 1);
optimal_input_vec = [repmat([-7;0.5], 10, 1);
                     repmat([15;15], 5, 1);
                     repmat([15;-0.5], 13, 1)];
time_horizon = length(optimal_input_vec)/2;

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
concat_state_realization = generateMonteCarloSims(n_monte_carlo, target_sys, ...
    target_init_state, time_horizon, optimal_input_vec);

% Target system support computation
concat_target_sys_dist_poly = ...
    Polyhedron('lb', repmat([dist_min;dist_min], time_horizon, 1), ...
               'ub', repmat([dist_max;dist_max], time_horizon, 1));
[Z, H, G] = target_sys.getConcatMats(time_horizon);
for tindx = 1:plot_t_skip:time_horizon-1
    fprintf('Plotting time: %d\n', tindx);
    relv_indx = [4*tindx+1, 4*tindx+3];
    support_poly = Z(relv_indx,:) * target_init_state + ...
        H(relv_indx, :) * optimal_input_vec + ...
        G(relv_indx, :) * concat_target_sys_dist_poly;
    plot(support_poly, 'alpha', 0.2, 'color', 'y');
    scatter(concat_state_realization(relv_indx(1)+4,1:skip_monte_carlo:end), ...
        concat_state_realization(relv_indx(2)+4,1:skip_monte_carlo:end), ...
        'ro', 'filled');
    drawnow
end

% Target system distribution computation
% [set_of_polytopes, hpoly] = polytopesFromMonteCarloSims( ...
%      concat_state_realization, 4, [1,3], {'color', 'red','alpha',0.2}, ...
%      plot_t_skip);

% Plot mean trajectory
mean_trajectory = mean(concat_state_realization,2);
plot(mean_trajectory(1:4:end), mean_trajectory(3:4:end), 'r--', 'linewidth', 2);
                            

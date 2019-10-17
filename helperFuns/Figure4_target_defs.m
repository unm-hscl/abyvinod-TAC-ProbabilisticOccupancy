% Target dynamics
srtinit

target_init_state = [15;0;-0.5;0];
relv_states = [1,3];
n_monte_carlo = 1e3;
skip_mc = 1;                        %Skip MC particles
dist_peak = 0;
dist_delta = 1;
dist_min = dist_peak - dist_delta;
dist_max = dist_peak + dist_delta;
dist_mean = dist_peak;              % Because of symmetry
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
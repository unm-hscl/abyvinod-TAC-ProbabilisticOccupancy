% Target dynamics
srtinit

n_monte_carlo = 1e4;
skip_mc = 1e2;                        %Skip MC particles
plot_t_skip = 3;

% %% Target system definition
% % Double integrator-based approach
% relv_states = [1,3];
% dist_peak = 0;
% dist_delta = 1.5;
% dist_min = dist_peak - dist_delta;
% dist_max = dist_peak + dist_delta;
% dist_mean = dist_peak;              % Because of symmetry
% target_init_state = [15;0;-0.5;0];
% target_affine_vec = [repmat([-7;0.5], 10, 1);
%                      repmat([15;15], 5, 1);
%                      repmat([15;-0.5], 10, 1)];
% time_horizon = length(target_affine_vec)/2;
% sys_DI_1D = getChainOfIntegLtiSystem(2, sampling_time, Polyhedron('lb',-1,'ub',1));
% target_sys_state_mat = blkdiag(sys_DI_1D.state_mat, sys_DI_1D.state_mat);
% target_sys_input_mat = blkdiag(sys_DI_1D.input_mat, sys_DI_1D.input_mat);
% target_sys_dist_mat = blkdiag(sys_DI_1D.input_mat, sys_DI_1D.input_mat);
% pd = makedist('Triangular','a', dist_min,'b', dist_peak,'c', dist_max);
% target_sys_dist = RandomVector('UserDefined', @(N) random(pd,2,N));
% target_sys = LtiSystem('StateMatrix', target_sys_state_mat, ...
%                        'InputMatrix', target_sys_input_mat, ...
%                        'DisturbanceMatrix', target_sys_dist_mat, ...
%                        'Disturbance', target_sys_dist, ...
%                        'InputSpace', ...
%                         Polyhedron('lb', -[inf,inf], 'ub', [inf,inf]));
% concat_target_sys_dist_poly = ...
%     Polyhedron('lb', repmat([dist_min;dist_min], time_horizon, 1), ...
%                'ub', repmat([dist_max;dist_max], time_horizon, 1));
% target_concat_mean_dist_vec = repmat([dist_mean; dist_mean], time_horizon, 1);                         

% Dubins vehicle dynamics
target_relv_states = [1,2];
sampling_time = 0.05;                           % Sampling time
initial_heading = -pi;                       % Initial heading 
target_init_state = [15;-0.5];
turning_rate_seq = [0*ones(1,3), ...
                    -5*ones(1,3), ...
                    -10*ones(1,3), ...
                    -5*ones(1,2), ...
                    0*ones(1,14)];
time_horizon = length(turning_rate_seq);
dist_delta = 5;
dist_peak = 20;
dist_min = dist_peak - dist_delta;
dist_max = dist_peak + dist_delta;
dist_mean = dist_peak;
v_rv_pdf_obj = makedist('Triangular','a', dist_min,'b', dist_peak, ...
    'c', dist_max);
v_rv = RandomVector('UserDefined', @(N) v_rv_pdf_obj.random(1,N));
target_sys = getDubinsCarLtv('vel-dist', turning_rate_seq', initial_heading, ...
    sampling_time, v_rv); 
target_affine_vec = zeros(0,1);
concat_target_sys_dist_poly = ...
    Polyhedron('lb', dist_min * ones(time_horizon, 1), ...
               'ub', dist_max * ones(time_horizon, 1));
target_concat_mean_dist_vec = dist_mean * ones(time_horizon, 1);

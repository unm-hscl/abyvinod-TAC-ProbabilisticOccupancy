%% Use fmincon for constrained optimization for permitted intercept zones
% Three steps:
% 1. Compute feas_poly as the intersection of the following polytopes:
%    a. Pursuer_i forward reach set
%    b. Support of the target
%    c. Intercept allowed regions
% 2. Compute an initial solution to the problem via a convex optimization
% problem
%    Project the mean target position to the feas_poly.
% 3. Use fmincon to solve the problem from this starting point
%    You can choose to use the gradient or not
use_gradient = false;
fmincon_optimoptions = optimoptions('fmincon', 'display', 'iter', ...
    'SpecifyObjectiveGradient', use_gradient);
max_feas_prob_catch_location = zeros(2, count_feas);
max_feas_prob_catch_value = zeros(1, count_feas);
init_feas_prob_catch_location = zeros(2, count_feas);
init_feas_prob_catch_value = zeros(1, count_feas);
elapsed_time_cvx = zeros(1, count_feas);
elapsed_time_fmincon = zeros(1, count_feas);

for feas_list_indx = 1:count_feas    
    % Get the parameters for feasible polytope
    pursuer_indx = feas_list(feas_list_indx, 1);
    t_indx_plus1 = feas_list(feas_list_indx, 2);   % Time goes from 0 to time_horizon
    pursuer_cvx_indx = feas_list(feas_list_indx, 3);
    fprintf('%2d/%2d. Computing for time %d\n', feas_list_indx, count_feas, ...
        t_indx_plus1 - 1);
    
    % Get the feasible polytope
    % We use t_indx_plus1 directly since t=0 sets are part of these lists
    feas_poly = pursuer_interceptable_position_set(pursuer_indx, ...
            t_indx_plus1, pursuer_cvx_indx).intersect( ...
                target_support_position(t_indx_plus1));
    feas_poly_outerapprox = feas_poly.outerApprox;
    
    % Compute target mean position (1 maps to 0) 
    % Z, H, G does not have initial state
    % -2 because t_indx_plus1 -1 is the correct time and an additional -1
    % ensures that relv_states (added on top of it) is of the correct time
    % snapshot
    if t_indx_plus1 > 1
        relv_indx = 4*(t_indx_plus1-2) + relv_states;
        target_mean_position = target_Z(relv_indx,:)*target_init_state+ ...
            target_H(relv_indx,:)*target_affine_vec + ...
            target_G(relv_indx,:)*repmat([dist_mean; dist_mean], time_horizon, 1);
    else
        target_mean_position = target_init_state;
    end
    
    
    % Compute initial state via projection (CVX)
    timer = tic;
    cvx_begin quiet
        variable initial_guess(2,1);
        
        minimize (norm(initial_guess - target_mean_position))
        subject to
            feas_poly.A * initial_guess <= feas_poly.b;
    cvx_end
    elapsed_time_cvx(feas_list_indx) = toc(timer);
    
    init_feas_prob_catch_location(:, feas_list_indx) = initial_guess;
    % t_indx_plus1 -1 because t_indx_plus1 starts from 0 and here the
    % correct t_indx is required
    target_affine_vec_slice = target_affine_vec(1:2*(t_indx_plus1-1));
    catch_prob = @(x) Figure4_occupy_fun_Levi(x + catch_box, t_indx_plus1-1, ...
        target_sys, relv_states, target_init_state, target_affine_vec_slice, ...
        dist_delta, dist_peak);

    timer = tic;
    init_feas_prob_catch_value(feas_list_indx) = catch_prob(initial_guess);
    if init_feas_prob_catch_value(feas_list_indx) >= 1 - zero_catch_prob
        max_feas_prob_catch_location(:, feas_list_indx) = initial_guess;
        max_feas_prob_catch_value(feas_list_indx) = ...
            init_feas_prob_catch_value(feas_list_indx);
        disp('Initial guess returned a probability of one');
    else
        fprintf('Initial guess returned a probability of %1.4f\n', ...
            init_feas_prob_catch_value(feas_list_indx));
        % Define the objective for fmincon
        if ~use_gradient
            % t_indx_plus1 -1 because t_indx_plus1 starts from 0 and here the
            % correct t_indx is required
            obj = @(x) -log(Figure4_occupy_fun_Levi(...
                x + catch_box, t_indx_plus1 - 1, target_sys, relv_states, ...
                target_init_state, target_affine_vec_slice, dist_delta, ...
                dist_peak, false));
        else
            % t_indx_plus1 -1 because t_indx_plus1 starts from 0 and here the
            % correct t_indx is required
            obj = @(x) obj_with_grad(x, catch_box, t_indx_plus1 - 1, target_sys, ...
                relv_states, target_init_state, target_affine_vec_slice, ...
                dist_delta, dist_peak);
        end
        % Compute the best capture location            
        [max_feas_prob_catch_location(:, feas_list_indx), F] = fmincon(...
            obj, initial_guess, feas_poly.A, feas_poly.b, [],...
            [], [], [], [], fmincon_optimoptions);
        fprintf('Catch probability: %1.4f\n', exp(-F));
        max_feas_prob_catch_value(feas_list_indx) = exp(-F);         
    end 
    elapsed_time_fmincon(feas_list_indx) = toc(timer);
    fprintf('Elapsed time: %1.2f s (cvx+fmincon)\n', ...
        elapsed_time_cvx(feas_list_indx) + elapsed_time_fmincon(feas_list_indx));
end

% Requires pursuer_sys, time_horizon, relv_states
% Computes position sets (projection of reach sets to x,y)
% Time goes from 0 to time_horizon

pursuer_mat_filename = 'pursuit_pursuer_position_set.mat';
elapsed_time_pursuer_reach = zeros(time_horizon + 1, 1);

if exist(pursuer_mat_filename, 'file')
    fprintf('### LOADED position sets from %s\n', pursuer_mat_filename);
    load(pursuer_mat_filename);
else
    pursuer_position_sets_zero_state_unit_input = ones(2,0)*Polyhedron();
    % Equal limits on x,y acceleration inputs
    pursuer_concat_lower_limit = -repmat([1;1], time_horizon, 1);
    pursuer_concat_upper_limit = -pursuer_concat_lower_limit;
    pursuer_concat_unit_input_space = ...
        Polyhedron('lb', pursuer_concat_lower_limit, ...
                   'ub', pursuer_concat_upper_limit);
    [pursuer_Z, pursuer_H, pursuer_G] = pursuer_sys.getConcatMats(time_horizon);
    for t_indx = 1:time_horizon
        fprintf('Computing pursuer forward reach (position) set for t=%d\n', ...
            t_indx);
        timer = tic;
        temp_poly = pursuer_H(4*(t_indx-1) + relv_states,:) * ...
            pursuer_concat_unit_input_space;
        temp_poly.minHRep();
        % t_indx + 1 to accommodate space for t=0
        elapsed_time_pursuer_reach(t_indx + 1) = toc(timer);
        pursuer_position_sets_zero_state_unit_input(t_indx) = temp_poly;
    end
    % Add the t=0 case
    timer = tic;
    pursuer_position_sets_zero_state_unit_input = ...
        [ones(2,0)*Polyhedron(), pursuer_position_sets_zero_state_unit_input];
    elapsed_time_pursuer_reach(1) = toc(timer);    
    save(strcat('./helperFuns/',pursuer_mat_filename), 'pursuer_Z', ...
        'pursuer_position_sets_zero_state_unit_input', ...
        'elapsed_time_pursuer_reach');
    fprintf('### SAVED position sets in %s\n', pursuer_mat_filename);
end
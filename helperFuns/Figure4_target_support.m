% Requires target_sys, time_horizon, relv_states, dist_min, dist_max,
% target_init_state, target_affine_vec
%
% Computes position sets (projection of reach sets to x,y)
% Time goes from 0 to time_horizon

target_mat_filename = 'Figure4_target_support.mat';

if exist(target_mat_filename, 'file')
    fprintf('### LOADED position supports from %s\n', target_mat_filename);
    load(target_mat_filename);
else
    concat_target_sys_dist_poly = ...
    Polyhedron('lb', repmat([dist_min;dist_min], time_horizon, 1), ...
               'ub', repmat([dist_max;dist_max], time_horizon, 1));
    [target_Z, target_H, target_G] = target_sys.getConcatMats(time_horizon);
    target_support_position = [ones(2,0) * Polyhedron()];
    for t_indx = 1:time_horizon
        fprintf('Computing support (position) for time: %d\n', t_indx);        
        relv_indx = 4*(t_indx-1) + relv_states;
        target_support_position(t_indx) = ...
            target_Z(relv_indx,:) * target_init_state + ...
            target_H(relv_indx, :) * target_affine_vec + ...
            target_G(relv_indx, :) * concat_target_sys_dist_poly;
    end
    target_support_position = [ones(2,0) * Polyhedron(), target_support_position];

    save(strcat('./helperFuns/',target_mat_filename), ...
        'target_support_position', 'target_Z', 'target_H', ...
        'target_G');
    fprintf('### SAVED position supports in %s\n', target_mat_filename);
end
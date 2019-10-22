% Requires
% time_horizon
% relv_states

pursuer_relv_states = [1,3];
pursuer1_initial_state = [8;0;0.5;0];
pursuer2_initial_state = [15;0;10;0];
pursuer3_initial_state = [27;0;1.5;0];
pursuer_u_limit = 4;

% Pursuer system definitions
sys_DI_1D = getChainOfIntegLtiSystem(2, sampling_time, ...
    Polyhedron('lb',-1,'ub',1));
pursuer_sys_state_mat = blkdiag(sys_DI_1D.state_mat, sys_DI_1D.state_mat);
pursuer_sys_input_mat = blkdiag(sys_DI_1D.input_mat, sys_DI_1D.input_mat);
pursuer_input_space = Polyhedron('lb',-[1;1],'ub',[1;1]);
pursuer_sys = LtiSystem('StateMatrix', pursuer_sys_state_mat, ...
                        'InputMatrix', pursuer_sys_input_mat, ...
                        'InputSpace', pursuer_input_space);

% Obtain the pursuer reach set with zero state and unit input                        
pursuit_pursuer_unit_input_position_set

% elapsed_time_pursuer_reach
if ~exist('elapsed_time_pursuer_reach','var')
    throw('Expected elapsed_time_pursuer_reach to be defined from '+ ...
        'Figure4_pursuer_unit_input_position_set');
end

% Compute the zero input (natural dynamics) of the pursuers
pursuer_position_set_1_zero_input = zeros(2,time_horizon);
pursuer_position_set_2_zero_input = zeros(2,time_horizon);
pursuer_position_set_3_zero_input = zeros(2,time_horizon);
for t_indx=1:time_horizon
    pursuer_position_set_1_zero_input(:, t_indx) = ...
        pursuer_Z(pursuer_sys.state_dim * (t_indx-1) + pursuer_relv_states,:)...
            * pursuer1_initial_state;
    pursuer_position_set_2_zero_input(:, t_indx) = ...
        pursuer_Z(pursuer_sys.state_dim * (t_indx-1) + pursuer_relv_states,:)...
            * pursuer2_initial_state;
    pursuer_position_set_3_zero_input(:, t_indx) = ...
        pursuer_Z(pursuer_sys.state_dim * (t_indx-1) + pursuer_relv_states,:)...
            * pursuer3_initial_state;
end

% Plot the forward position sets for the pursuers
pursuer_team_position_set_zero_input = [
    pursuer1_initial_state(pursuer_relv_states,1), pursuer_position_set_1_zero_input;
    pursuer2_initial_state(pursuer_relv_states,1), pursuer_position_set_2_zero_input;
    pursuer3_initial_state(pursuer_relv_states,1), pursuer_position_set_3_zero_input];
pursuer_interceptable_position_set = [ones(2,0)*Polyhedron()];
for t_indx_plus1 = 1:time_horizon+1
    % Time goes from 0 to time_horizon for both
    %   pursuer_position_set_zero_input and 
    %   pursuer_position_sets_zero_state_unit_input
    pursuer_reach_timer = tic;
    for pursuer_indx = 1:3
        % A 2x(time_horizon + 1) matrix of pursuer positions (under natural
        % dynamics)
        pursuer_position_set_zero_input = ...
            pursuer_team_position_set_zero_input(...
                2*(pursuer_indx-1)+1 : 2*(pursuer_indx-1)+2, :);
        for poly_indx = 1:3
            % Get the forward reach set
            if ~pursuer_position_sets_zero_state_unit_input(t_indx_plus1).isEmptySet()
                pursuer_actual_position_set = ...
                    pursuer_position_set_zero_input(:, t_indx_plus1) + ...
                        pursuer_u_limit * ...
                        pursuer_position_sets_zero_state_unit_input(t_indx_plus1);
            else
                pursuer_actual_position_set = Polyhedron('V', ...
                    pursuer_position_set_zero_input(:, t_indx_plus1)');
            end
            % Forward reach set at time t intersected with convex intercept zone
            pursuer_interceptable_position_set( ...
                            pursuer_indx, t_indx_plus1, poly_indx) = ...
                pursuer_actual_position_set.intersect( ...
                    pursuer_cvx(pursuer_indx, poly_indx));
        end
    end
    elapsed_time_pursuer_reach(t_indx_plus1) = elapsed_time_pursuer_reach(t_indx_plus1) +... 
        toc(pursuer_reach_timer);
end
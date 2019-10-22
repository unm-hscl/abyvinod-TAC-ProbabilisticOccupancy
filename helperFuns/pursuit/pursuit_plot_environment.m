% Plot
if plot_layout_only
    figure(1);clf;
else
    figure(4);clf;
end
plot(area_of_interest, 'color', 'b', 'alpha',0.2);
hold on;
plot(pursuer_cvx(1,1), 'color', 'g', 'alpha',0.4);
plot(pursuer_cvx(1,2), 'color', 'g', 'alpha',0.4);
plot(pursuer_cvx(1,3), 'color', 'g', 'alpha',0.4);
plot(pursuer_cvx(2,1), 'color', 'g', 'alpha',0.4);
plot(pursuer_cvx(2,2), 'color', 'g', 'alpha',0.4);
plot(pursuer_cvx(2,3), 'color', 'g', 'alpha',0.4);
plot(pursuer_cvx(3,1), 'color', 'g', 'alpha',0.4);
plot(pursuer_cvx(3,2), 'color', 'g', 'alpha',0.4);
plot(pursuer_cvx(3,3), 'color', 'g', 'alpha',0.4);
axis equal;
box on;
set(gca,'FontSize',fontSize*2);
grid on;
xlim([-2,37]);
ylim([-2,15]);
xlabel('x','Interpreter','latex');
ylabel('y','Interpreter','latex');


if plot_layout_only
    for t_indx_plus1 = 2:plot_t_skip:time_horizon+1
        % Time goes from 0 to time_horizon for both
        %   target_support_position and target_concat_state_realization
        fprintf('Plotting time support: %d\n', t_indx_plus1-1);
        plot(target_support_position(t_indx_plus1), 'alpha', 0.2, 'color', 'y');
        relv_indx = target_sys.state_dim*(t_indx_plus1-1) + target_relv_states;
        scatter(target_concat_state_realization(relv_indx(1),1:skip_mc:end), ...
                target_concat_state_realization(relv_indx(2),1:skip_mc:end), ...
                'ro', 'filled');
        drawnow
    end
    for pursuer_indx = 1:3
        pursuer_initial_position = pursuer_team_position_set_zero_input( ...
            2*(pursuer_indx-1) + [1,2],1);
        scatter(pursuer_initial_position(1), pursuer_initial_position(2), ...
                1,'bx', 'LineWidth', 5);        
    end
    % Plot mean trajectory
    mean_trajectory = mean(target_concat_state_realization,2);
    plot(mean_trajectory(target_relv_states(1):target_sys.state_dim:end), ...
         mean_trajectory(target_relv_states(2):target_sys.state_dim:end), 'r--', ...
         'linewidth', 2);
else
    for t_indx_plus1=time_horizon+1:-1:1
        for pursuer_indx = 1:3
            for poly_indx = 1:3
                if ~pursuer_interceptable_position_set( ...
                                pursuer_indx, t_indx_plus1, poly_indx).isEmptySet() && ...
                        abs(mod(t_indx_plus1, plot_t_skip))<1e-8
                    fprintf('Plotting time for pursuer %d: %d\n', pursuer_indx, ...
                        t_indx_plus1-1);
                    plot(pursuer_interceptable_position_set( ...
                                pursuer_indx, t_indx_plus1, poly_indx), ...
                        'alpha', 0.6, 'color', 'c');
                end
            end
        end
    end
    % Pursuer best intercept location
    t_plus1_vec = [];
    for pursuer_indx = 1:3
        pursuer_initial_position = pursuer_team_position_set_zero_input( ...
            2*(pursuer_indx-1) + [1,2],1);
        pursuer_valid_indx = find(feas_list(:,1) == pursuer_indx);
        prob_values = max_feas_prob_catch_value(pursuer_valid_indx);
        [max_prob, prob_indx] = max(prob_values);
        best_feas_list_arg = pursuer_valid_indx(prob_indx);
        
        opt_intercept_t = feas_list(best_feas_list_arg,2) - 1;
        t_plus1_vec = [t_plus1_vec, opt_intercept_t + 1];
        
        [pursuer_Z, pursuer_H, ~] = pursuer_sys.getConcatMats(opt_intercept_t);
        % Compute the controller
        cvx_begin quiet
            variable U(pursuer_sys.input_dim * opt_intercept_t, 1)
            variable X(pursuer_sys.state_dim * opt_intercept_t, 1)
            
            minimize (norm(U) + norm(X(end - pursuer_sys.state_dim+1:end)))
            
            subject to
                X == pursuer_Z * [pursuer_initial_position(1);0;
                    pursuer_initial_position(2);0] + pursuer_H * U;
                X(end - pursuer_sys.state_dim + pursuer_relv_states(1)) == ...
                    max_feas_prob_catch_location(1, best_feas_list_arg);
                X(end - pursuer_sys.state_dim + pursuer_relv_states(2)) == ...
                    max_feas_prob_catch_location(2, best_feas_list_arg);
                abs(U) <= pursuer_u_limit;
        cvx_end
        pursuer_trajectory = reshape(X, pursuer_sys.state_dim, opt_intercept_t);
        target_relv_indx = target_sys.state_dim * opt_intercept_t + ...
            target_relv_states;
        scatter(pursuer_initial_position(1), pursuer_initial_position(2), ...
                200,'bx', 'LineWidth', 5);
        scatter(target_concat_state_realization(target_relv_indx(1),1:skip_mc:end), ...
                target_concat_state_realization(target_relv_indx(2),1:skip_mc:end), ...
                'ro', 'filled');
        scatter(pursuer_trajectory(pursuer_relv_states(1),:), ...
                pursuer_trajectory(pursuer_relv_states(2),:), 100, 'bo'); 
        scatter(max_feas_prob_catch_location(1, best_feas_list_arg), ...
                max_feas_prob_catch_location(2, best_feas_list_arg), ...
                200, 'bo', 'filled');
        pursuer_catch_box = catch_box + ...
            max_feas_prob_catch_location(:, best_feas_list_arg);
        plot(pursuer_catch_box, 'alpha', 0.3, 'color', 'b')
        mc_prob_estim = sum(pursuer_catch_box.contains( ...
            target_concat_state_realization(target_relv_indx,1:skip_mc:end)))...
            /n_monte_carlo;
        fprintf('Prob: %1.4f | MC Estim. Prob: %1.4f\n', max_prob, mc_prob_estim);
    end    
end
clearvars;
clc;
close all;

addpath('./helperFuns/Fig789_pursuit/');
%% Global parameters
fontSize = 30;
plot_layout_only = 1;

%% Construct the approximate map
pursuit_polytope_defs

%% Construct the target, forward reach sets, and distribution
pursuit_target_defs

%% Test: Compute probability
% query_point = [19.5;10.5];
% my_eps = 2e0;
% time_step = 25;
% unit_box = my_eps * Polyhedron('lb',-[1;1],'ub',[1;1]);
% query_box = query_point + unit_box;
% plot(query_box, 'color','b', 'alpha',0.3);
% prob = pursuit_occupy_fun_Levi(query_box, time_step, target_sys, [1,3], ...
%         target_init_state, optimal_input_vec(1:2*time_step), dist_delta, ...
%         dist_peak);
% % Monte-Carlo-simulation based validation    
% relv_sims = concat_state_realization(target_sys.state_dim * time_step + ...
%     target_relv_states,:);    
% prob_mcarlo = sum(query_box.contains(relv_sims))/n_monte_carlo;
% fprintf('Probability : %1.4f | MonteCarlo probability : %1.4f\n', ...
%     prob, prob_mcarlo);    

%% Construct the pursuer forward reach sets
pursuit_pursuer_defs

%% Target system support computation
% Get the target support
pursuit_target_support

% Generate Monte Carlo simulation                    
target_concat_state_realization = generateMonteCarloSims(n_monte_carlo, ...
    target_sys, target_init_state, time_horizon, target_affine_vec);

%% Plot environment
pursuit_plot_environment;

% Stop here if you want to see the layout 
%% Define catch probability
catch_box_half_length = 5e-1;
zero_catch_prob = 1e-4;
catch_box = catch_box_half_length * Polyhedron('lb',-[1;1],'ub',[1;1]);

%% Find non-empty feasible intersect locations
feasible_intercept_locations = [ones(2,0) * Polyhedron()];
count_infeas = 0;
count_feas = 0;
feas_list = [];       % pursuer_indx, t_indx_plus1, pursuer_cvx_indx
elapsed_time_feas_poly = zeros(time_horizon + 1, 1);
feas_list_polytope = [];
for t_indx_plus1 = 1:time_horizon+1
    % Time goes from 0 to time_horizon for both
    %   pursuer_position_set_zero_input and 
    %   pursuer_position_sets_zero_state_unit_input
    target_support_poly_plus_box = target_support_position(t_indx_plus1) + ...;
        catch_box;
    feas_poly_timer=tic;
    for pursuer_indx = 1:3
        for pursuer_cvx_indx = 1:3
            temp_poly = pursuer_interceptable_position_set(pursuer_indx, ...
                t_indx_plus1, pursuer_cvx_indx);
            temp_poly = temp_poly.intersect(target_support_poly_plus_box);
            if temp_poly.isEmptySet()
                count_infeas = count_infeas + 1;
            else
                count_feas = count_feas + 1;
                feas_list = [feas_list;pursuer_indx, t_indx_plus1, pursuer_cvx_indx];
                feas_list_polytope = [feas_list_polytope, temp_poly];
            end                
        end
    end
    elapsed_time_feas_poly(t_indx_plus1) = toc(feas_poly_timer);
end
fprintf('Need to solve the catch problem with %d polytopes (%1.3f %%)\n', ...
    count_feas, count_feas/(3*3*time_horizon)*100);

%% Use fmincon for constrained optimization for permitted intercept zones
pursuit_optimize_via_fmincon

%% Plot the mean times
pursuit_plot_results

%% Final touches to the plot
plot_layout_only = 0;
pursuit_plot_environment;
figure(1);
axis equal;
xlim([-2,37]);
ylim([-2,15]);
% rmpath('./helperFuns/Fig789_pursuit/');

%% Save figures
figure(1);
saveas(gcf, 'figs/PursuitLayout.png', 'png');
savefig(gcf, 'figs/PursuitLayout.fig', 'compact');
figure(4);
saveas(gcf, 'figs/PursuitSolve.png', 'png');
savefig(gcf, 'figs/PursuitSolve.fig', 'compact');
figure(3);
saveas(gcf, 'figs/PursuitCatchProb.png', 'png');
savefig(gcf, 'figs/PursuitCatchProb.fig', 'compact');


%% Results
fprintf('Solve time for %d optimization problems (fmincon + cvx): % 1.2f (%1.1f minutes)\n', ...
    length(elapsed_time_cvx), sum(elapsed_time_fmincon + elapsed_time_cvx), ...
    sum(elapsed_time_fmincon + elapsed_time_cvx)/60);
fprintf('Total computation time: % 1.2f (%1.1f minutes)\n', sum(total_time), ...
    sum(total_time)/60);

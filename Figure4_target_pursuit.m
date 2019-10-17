clearvars;
clc;

addpath('./helperFuns/');

%% Construct the approximate map
Figure4_polytope_defs

%% Construct the target, forward reach sets, and distribution
Figure4_target_defs

%% Test: Compute probability
% query_point = [19.5;10.5];
% my_eps = 2e0;
% time_step = 25;
% unit_box = my_eps * Polyhedron('lb',-[1;1],'ub',[1;1]);
% query_box = query_point + unit_box;
% plot(query_box, 'color','b', 'alpha',0.3);
% prob = Figure4_occupy_fun_Levi(query_box, time_step, target_sys, [1,3], ...
%         target_init_state, optimal_input_vec(1:2*time_step), dist_delta, ...
%         dist_peak);
% % Monte-Carlo-simulation based validation    
% relv_sims = concat_state_realization(4*time_step + [1;3],:);    
% prob_mcarlo = sum(query_box.contains(relv_sims))/n_monte_carlo;
% fprintf('Probability : %1.4f | MonteCarlo probability : %1.4f\n', ...
%     prob, prob_mcarlo);    

% %% Construct the pursuer forward reach sets
% Figure4_pursuer_defs

%% Construct time-stamped reach sets for each of the pursuers

%% Use fmincon for constrained optimization for permitted intercept zones

figure(1);
axis equal;
xlim([-2,37]);
ylim([-3,19]);

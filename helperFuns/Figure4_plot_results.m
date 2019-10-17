elapsed_time_solver = elapsed_time_fmincon;

% Decode the feas_list into indexed time stamps
enum_feas_list = [feas_list, (1:count_feas)'];

pursuer_compute_time_cvx = zeros(time_horizon+1,1);
pursuer_compute_time_solver = zeros(time_horizon+1,1);
% pursuer_compute_time_feas_check;

for pursuer_indx = 1:3
    pursuer_indx_match = enum_feas_list(:,1) == pursuer_indx;
    time_and_indx = [enum_feas_list(pursuer_indx_match,2), ...
                     enum_feas_list(pursuer_indx_match,4)];
    unique_time_steps = unique(time_and_indx(:,1))';
    for t_indx_plus1 = unique_time_steps
        indices_to_use = time_and_indx(time_and_indx(:,1) == t_indx_plus1,2);
        
        % Pursuer cvx time
        pursuer_compute_time_cvx(t_indx_plus1) = ...
            pursuer_compute_time_cvx(t_indx_plus1) + ...
            sum(elapsed_time_cvx(indices_to_use));

        % Pursuer solver time
        pursuer_compute_time_solver(t_indx_plus1) = ...
            pursuer_compute_time_solver(t_indx_plus1) + ...
            sum(elapsed_time_solver(indices_to_use));

        if length(indices_to_use) > 1
            fprintf('Pursuer %d | Multiple indices (%d) at t=%d | %s \n', ...
                pursuer_indx, length(indices_to_use), t_indx_plus1-1, ...
                num2str(indices_to_use'));            
            
        else
            fprintf('Pursuer %d | Single index at t=%d         | %d \n', ...
                pursuer_indx, t_indx_plus1-1, indices_to_use);
        end        
    end
end


%% Plot stem plot of computational time
figure(2);
clf;
hold on;
computeTime_markerSize = 12.5;
computeTime_linewidth = 3;
% stem(0:time_horizon, elapsed_time_pursuer_reach, ...
%     'DisplayName', 'Pursuer reach set');
% stem(0:time_horizon, elapsed_time_target_support, ...
%     'DisplayName', 'Target support');
% stem(0:time_horizon, elapsed_time_feas_poly, ...
%     'DisplayName', 'Feasibility check');
stem(0:time_horizon, elapsed_time_feas_poly + elapsed_time_pursuer_reach + ...
    elapsed_time_target_support, 'ms', 'filled', ...
    'MarkerSize', computeTime_markerSize, 'LineWidth', computeTime_linewidth,...
    'DisplayName', 'Feasibility check');
stem(0:time_horizon, pursuer_compute_time_cvx, 'bo', 'filled', ...
    'MarkerSize', computeTime_markerSize, 'LineWidth', computeTime_linewidth,...
    'DisplayName', 'Initial guess');
stem(0:time_horizon, pursuer_compute_time_solver, 'rv', 'filled', ...
    'MarkerSize', computeTime_markerSize, 'LineWidth', computeTime_linewidth,...
    'DisplayName', 'Optimization');
total_time = elapsed_time_feas_poly + elapsed_time_pursuer_reach + ...
    elapsed_time_target_support + pursuer_compute_time_cvx + ...
    pursuer_compute_time_solver;
stem(0:time_horizon, total_time, 'kp', 'filled',  ...
    'MarkerSize', computeTime_markerSize, 'LineWidth', computeTime_linewidth,...
    'DisplayName', 'Total time');
leg = legend();
set(leg, 'Location','Best');
grid on;
set(gca,'FontSize', fontSize);
box on;

% Plot reach probability as a side-by-side plot
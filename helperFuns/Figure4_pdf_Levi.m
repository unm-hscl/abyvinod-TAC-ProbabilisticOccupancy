function pdf = Figure4_pdf_Levi(query_point, half_box_side, time_step, ...
        target_sys, relv_states, target_init_state, optimal_input_vec, ...
        dist_delta, dist_peak)
    % First argument is the point at which pdf is sought
    % Second argument is the half of the box side length used for
    % Lebesgue's density theorem
    % Other arguments coincide with Figure4_occupy_fun_Levi
    query_box = half_box_side * Polyhedron('lb',[-1;-1],'ub',[1;1]) + ...
        query_point;
    prob = Figure4_occupy_fun_Levi(query_box, time_step, target_sys, ...
        relv_states, target_init_state, optimal_input_vec(1:2*time_step), ...
        dist_delta, dist_peak);
    query_box_volume = (2 * half_box_side)^2;
    pdf = prob/query_box_volume;
end
function prob = Figure4_occupy_fun_Levi(query_box, time_step, sys, ...
        relv_states, initial_state, optimal_input_vec, dist_delta, dist_peak)
    % query_box         - Probability of lying within this box is computed
    % time_step         - Time step into the future
    % relv_states       - components of the state space interested
    % initial_state     - Initial state of the system
    % optimal_input_vec - Affine term that will be multiplied with H
    % dist_delta        - Disturbance width (b-a in Wikipedia)
    % dist_peak         - Disturbance peak (c in Wikipedia)
    
    % Throw a warning if probability outside [-my_zero, 1+my_zero]
    my_zero = 1e-6;    
    bounds_on_integral = 400;
    [Z, H, G] = sys.getConcatMats(time_step); 
    
    % x_t = G W + (Z*x_0 + H * U) last 4 rows = d + T*W
    state_dim = sys.state_dim;
    d = Z(end-state_dim+relv_states,:) * initial_state + ...
        H(end-state_dim+relv_states,:) * optimal_input_vec;
    T = G(end-state_dim+relv_states,:);
    
    
    % Extract corners from query_box
    lb = min(query_box.V)';
    ub = max(query_box.V)';     
    my_integrand = @(t1, t2) integrand(t1, t2, d, T, dist_delta, dist_peak, ...
        lb, ub);
    % Do the integral
    prob_t2_integrated= @(t1) integral(@(t2) my_integrand(t1,t2), ...
        -bounds_on_integral, bounds_on_integral);
    prob = integral(@(t1) prob_t2_integrated(t1), -bounds_on_integral, ...
        bounds_on_integral, 'ArrayValued', true);
    
    if abs(imag(prob)) >= 1e-8
        throw('Got a complex-valued prob!')
    else
        prob = real(prob);
    end
    
    if prob < - my_zero || prob > 1 + my_zero
        warning('Got probability outside the range');
        keyboard
    else
        % Sanitize
        prob = max(min(prob, 1), 0);
    end
end

function val = integrand(t1, t2, d, T, dist_delta, dist_peak, lb, ub)
    t1 = repmat(t1, size(t2,1), size(t2,2));
    
    % Define the multiplication term
    val = zeros(size(t1));
    for indx_x = 1:size(t1,1)
        for indx_y = 1:size(t1,2)
            t_vec = [t1(indx_x, indx_y);
                     t2(indx_x, indx_y)];
            % Define the integrand 
            val(indx_x, indx_y) = extra_terms_levi(t_vec, lb, ub) * ...
                cfun_W(T'*t_vec, dist_delta, dist_peak) * ...
                exp(1i * d' * t_vec);
        end
    end
    if any(any(isnan(val))) || any(any(isinf(val)))
        fprintf('IsNan: %d | IsInf: %d\n', any(any(isnan(val))), ...
            any(any(isinf(val))));
        keyboard
    end    
end

function val = extra_terms_levi(t_vec, lb, ub)
    % Expects t_vec in a 2-dimensional column
    box_center = (lb+ub)/2;
    delta_bounds = (ub - lb)/2;
    % t_vec is 2-dimensional column vector
    % 2 because exp(1i * delta_bound * t(i)) - exp(01i * delta_bound * t(i))
    val = exp(-1i * (box_center'*t_vec)) * ...
        prod(2 * (mysinc(delta_bounds.*t_vec) .* delta_bounds));    
    % Constant of integration in Levi's formula
    val = val/(2*pi)/(2*pi);
end

function val = cfun_W(t_vec, dist_delta, dist_peak)
    % Expects t_vec in a 2T-dimensional column
    % Expects dist_peak as a scalar
    % Expects dist_delta as a 2T-dimensional column or a scalar
    indiv_args = (t_vec .* dist_delta)/2;
    % exp( 1i * dist_peak * sum(t_vec)) denotes the product term
    val = prod(mysinc(indiv_args))^2 * exp(1i * dist_peak * sum(t_vec));
end

function val = mysinc(x)
    my_zero = 1e-6;
    val_with_NaN = sin(x)./x;
    val_with_NaN(abs(x) <= my_zero) = 1;
    val = val_with_NaN;
end
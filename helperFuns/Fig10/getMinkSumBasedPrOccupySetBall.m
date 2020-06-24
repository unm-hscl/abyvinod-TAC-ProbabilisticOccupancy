function polytopic_overapproximation = getMinkSumBasedPrOccupySetBall(...
        obstacle_mu, obstacle_sigma, obstacle_radius, alpha_value, N_des)

    % Sample equally spaced unit vectors
    Amatrix = spreadPointsOnUnitSphere(2, N_des, 0)';
    
    obstacle_volume = pi * obstacle_radius^2;
    probability_threshold = alpha_value/obstacle_volume;
        
    if det(obstacle_sigma) < 1e-8        
        [V, D] = eig(obstacle_sigma);
        support_1D = V(:, 2);
        std_dev = sqrt(D(2, 2));
        bound_1D = std_dev * ...
            sqrt(-2 * log(sqrt(2 * pi) * std_dev * probability_threshold));
        level_set_polytope_V = [obstacle_mu - bound_1D * support_1D, ...
                                obstacle_mu + bound_1D * support_1D];
        support_vector_evaluations_level_set = ...
            max(Amatrix * level_set_polytope_V, [], 2);
    else
        % Compute the ellipsoid level set (from CDC 2017 paper)
        % x^t (R^2)\Sigma^{-1}<=1 is the ellipsoid
        R_squared =-2*log(probability_threshold*sqrt(det(2*pi*obstacle_sigma)));
        ellipsoid_shape_matrix = obstacle_sigma * R_squared;
        
        support_ellipsoid_shape_matrix_sqrt = ...
            sqrt(diag(Amatrix * ellipsoid_shape_matrix * Amatrix'));
        support_vector_evaluations_level_set = Amatrix * obstacle_mu ...
            + support_ellipsoid_shape_matrix_sqrt;        
    end
    bmatrix = support_vector_evaluations_level_set + obstacle_radius;
    % Overapproximation via polytopes
    polytopic_overapproximation = Polyhedron('H', [Amatrix  bmatrix]);
end

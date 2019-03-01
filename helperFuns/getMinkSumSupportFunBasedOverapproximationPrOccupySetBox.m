function polytopic_overapproximation = ...
             getMinkSumSupportFunBasedOverapproximationPrOccupySetBox(...
                                                    obstacle_mu,...
                                                    obstacle_sigma,...
                                                    obstacle_polytope,...
                                                    alpha_value,...
                                                    obstacle_volume,...
                                                    no_of_vectors)

    % Setup
    Amatrix = zeros(no_of_vectors,2);
    bmatrix = zeros(no_of_vectors,1);

    % Sample equally spaced unit vectors
    theta_vector = linspace(0, 2 * pi, no_of_vectors + 1);
    Amatrix = [cos(theta_vector(1:end-1));
               sin(theta_vector(1:end-1))]';

    % Compute the ellipsoid level set (from CDC 2017 paper)
    % x^t (R^2)\Sigma^{-1}<=1 is the ellipsoid
    probability_threshold = alpha_value/obstacle_volume;
    R_squared = - 2 *  log(sqrt(det(2*pi * obstacle_sigma)) * probability_threshold);
    ellipsoid_shape_matrix = obstacle_sigma * R_squared;
    bmatrix = Amatrix * obstacle_mu  + sqrt(diag(Amatrix * ellipsoid_shape_matrix * Amatrix')) + obstacle_polytope.support(Amatrix');
   
    % Overapproximation via polytopes
    polytopic_overapproximation = Polyhedron('H', [Amatrix  bmatrix]);
end

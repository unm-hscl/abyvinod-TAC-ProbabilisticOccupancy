function polytopic_overapproximation = ...
             getMinkSumSupportFunBasedOverapproximationPrOccupySet(...
                                                    obstacle_mu,...
                                                    obstacle_sigma,...
                                                    obstacle_radius,...
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
    if det(obstacle_sigma)>1e-8
        R_squared = - 2 *  log(sqrt(det(2*pi * obstacle_sigma)) * probability_threshold);
        ellipsoid_shape_matrix = obstacle_sigma * R_squared;
        bmatrix = Amatrix * obstacle_mu  + sqrt(diag(Amatrix * ellipsoid_shape_matrix * Amatrix')) + obstacle_radius;
        polytopic_overapproximation = Polyhedron('H', [Amatrix  bmatrix]);
    else
        % Specialized for 2D and a ball of obstacle_radius radius
        [V,D] = eig(obstacle_sigma);
        direction_of_sigma = V(:,2);
        % Since the diameter is the largest chord, the occupancy function has
        % the highest value when restricted along V(:,1) [perpendicular to
        % V(:,2)] when the point of interest overline{y} lies on V(:,2).
        % Furthermore, the highest value the occupancy function will ever have
        % is when \overline{y} is at \mu by the symmetry proposition.
        % Hence, replace m( mathcal{O}(overline{0})) with [-obstacle_radius,
        % obstacle_radius] and use the 1D FSRPD along V(:,2) with standard
        % deviation sqrt(non-zero eigenvalue of Sigma) and mean zero (So that
        % after lifting to V(:,1) and translation to obstacle_mu, we recover the
        % embedded 1D Gaussian in 2D space.
        R_squared = - 2 *  log(sqrt(2*pi * D(2,2)) * alpha_value/(2*obstacle_radius));
        one_extreme = sqrt(D(2,2) * R_squared);
        gaussian_level_set = Polyhedron('V',...
                                            [obstacle_mu' + direction_of_sigma' * one_extreme;
                                             obstacle_mu' - direction_of_sigma' * one_extreme]);
        obstacle_polytope =...
        Polyhedron('V',obstacle_radius*[cos(linspace(0,2*pi,100))',sin(linspace(0,2*pi,100))']);
        polytopic_overapproximation = gaussian_level_set + obstacle_polytope;
    end
    % Overapproximation via polytopes
end

function [polytopic_overapproximation,...
          boundary_point_of_interest,...
          Gaussian_set_polytope] = ...
             getMinkSumBasedOverapproximationPrOccupySet(...
                                                    obstacle_mu,...
                                                    obstacle_sigma,...
                                                    obstacle_polytope,...
                                                    obstacle_volume,...
                                                    sampling_circle_radius,...
                                                    alpha_value,...
                                                    no_of_points_on_the_circle)

    % Setup
    boundary_point_of_interest = zeros(2,no_of_points_on_the_circle);
    Amatrix = zeros(no_of_points_on_the_circle,2);
    bmatrix = zeros(no_of_points_on_the_circle,1);

    % Sample points on the sampling circle
    sampling_circle_center = obstacle_mu;
    theta_vector = linspace(0, 2 * pi, no_of_points_on_the_circle + 1);
    points_on_the_sampling_circle = sampling_circle_center + ...
               sampling_circle_radius * [cos(theta_vector(1:end-1));
                                         sin(theta_vector(1:end-1))];

    % Compute the ellipsoid level set (from CDC 2017 paper)
    % x^t (R^2)\Sigma^{-1}<=1 is the ellipsoid
    probability_threshold = alpha_value/volume(obstacle_polytope);
    R_squared = - 2 *  log(sqrt(det(2*pi * obstacle_sigma)) * probability_threshold);
    ellipsoid_shape_sqrtm = chol(inv(obstacle_sigma * R_squared));

    for point_index = 1 : no_of_points_on_the_circle
        point_of_interest = points_on_the_sampling_circle(:, point_index);

        %% Project the points on to ellipsoid and obtain the corresponding supporting
        %%hyperplane via its dual variable
        % minimize |x - x_0|_2
        % s.t.     x^\top inv(obstacle_sigma * R_squared) * x <= 1 
        % (norm form is better for cvx)
        cvx_begin quiet
            cvx_precision best
            variable projection_point(2)
            variable y(2)
            dual variable hyperplane_vector
            minimize norm(y)
            subject to
               norm(ellipsoid_shape_sqrtm * (projection_point - obstacle_mu)) <= 1
               hyperplane_vector : projection_point - point_of_interest == y
        cvx_end

        % Translate the optimal svalue of the fmincon into physical points
        boundary_point_of_interest(:,point_index) = projection_point;
        % The optimal problem satisfying Slater's condition implies the optimal
        % dual objective would be equal to the optimal primal objective value 
        % hyperplane_vector' * (projection_point - point_of_interest) = 
        % || projection_point - point_of_interest ||
        % We were after the computation of (projection_point - point_of_interest)
        hyperplane_vector_normalized =...
                                    (hyperplane_vector)/norm(hyperplane_vector);

        % The supporting hyperplane is given by: (P(x) - x_0)^T (x-x_0)<0
        Amatrix(point_index,:) = hyperplane_vector_normalized';
        bmatrix(point_index) = hyperplane_vector_normalized'*projection_point;
    end
    % Create the polytope using the hyperplanes
    Gaussian_set_polytope=Polyhedron('H', [Amatrix  bmatrix]);
    % Add these two together
    polytopic_overapproximation = Gaussian_set_polytope + obstacle_polytope;
end

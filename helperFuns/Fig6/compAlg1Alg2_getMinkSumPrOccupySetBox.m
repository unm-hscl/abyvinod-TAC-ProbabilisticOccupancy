function polytopic_overapproximation = compAlg1Alg2_getMinkSumPrOccupySetBox(...
    obstacle_mu, obstacle_sigma, obstacle_polytope, alpha_value, ...
    obstacle_volume, no_of_vectors)
    % Implement Algorithm 2
    
    % Sample equally spaced unit vectors
    theta_vector = linspace(0, 2 * pi, no_of_vectors + 1);
    Amatrix = [cos(theta_vector(1:end-1));
               sin(theta_vector(1:end-1))]';

    % Compute the ellipsoid level set (from CDC 2017 paper)
    % x^t (R^2)\Sigma^{-1}<=1 is the ellipsoid
    probability_threshold = alpha_value/obstacle_volume;
    R_squared = - 2 *  log(sqrt(det(2*pi * obstacle_sigma)) * probability_threshold);
    ellipsoid_shape_matrix = obstacle_sigma * R_squared;
    support_ellipsoid_shape_matrix_sqrt = sqrt(diag(Amatrix * ellipsoid_shape_matrix * Amatrix'));

%   % Via linprog
%     support_matrix = zeros(no_of_vectors, 1);
%     for indx = 1:no_of_vectors
%         [~,support_matrix(indx)] = linprog(Amatrix(indx,:), obstacle_polytope.A, obstacle_polytope.b);
%     end

%   % Via MPT
%     support_matrix = obstacle_polytope.support(Amatrix');

%   % Via support function
    support_matrix = max(Amatrix * obstacle_polytope.V',[],2);
    
    bmatrix = Amatrix * obstacle_mu  + support_ellipsoid_shape_matrix_sqrt + support_matrix;
   
    % Overapproximation via polytopes
    polytopic_overapproximation = Polyhedron('H', [Amatrix  bmatrix]);
end

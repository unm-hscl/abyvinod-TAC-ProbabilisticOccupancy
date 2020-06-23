function [polytopic_overapproximation, polytopic_underapproximation,...
          boundary_point_of_interest] = ...
                compAlg1Alg2_getTightPolytopicOverapproximationPrOccupySet(...
                    obstacle_mu, obstacle_sigma, obstacle_lower_bound,...
                    obstacle_upper_bound, bounding_ball_radius,...
                    alpha_value, no_of_points_on_the_bounding_ball, ...
                    fminconsettings)
    % Implement Algorithm 1
        
    % Setup
    log_alpha_value = log(alpha_value);
    boundary_point_of_interest = zeros(2,no_of_points_on_the_bounding_ball);
    Amatrix = zeros(no_of_points_on_the_bounding_ball,2);
    bmatrix = zeros(no_of_points_on_the_bounding_ball,1);

    % Nonlinear constraint for fmincon:  log(\alpha) - log(\phi(z)) \leq 0
    % Defined below
    nonlcon = @(x) PrOccupySetConstraint(x,...
                                         obstacle_mu,...
                                         obstacle_sigma,...
                                         obstacle_lower_bound,...
                                         obstacle_upper_bound,...
                                         log_alpha_value);

    % Sample points on the sampling circle
    sampling_circle_center = obstacle_mu;
    theta_vector = linspace(0, 2 * pi, no_of_points_on_the_bounding_ball + 1);
    points_on_the_bounding_ball = sampling_circle_center + ...
               bounding_ball_radius * [cos(theta_vector(1:end-1));
                                       sin(theta_vector(1:end-1))];
   
    % Iterate over all points on the sampling circle
    for point_index = 1 : no_of_points_on_the_bounding_ball
        point_of_interest = points_on_the_bounding_ball(:, point_index);

        % Compute the projection via fmincon
        lb_fmincon = sampling_circle_center-bounding_ball_radius*ones(2,1);
        ub_fmincon = sampling_circle_center+bounding_ball_radius*ones(2,1);
        [projection_point, ~] = fmincon(@(x) norm(x - point_of_interest),...
                                        obstacle_mu, [],[],[],[],...
                                        lb_fmincon, ub_fmincon,...
                                        nonlcon, fminconsettings);
        % Translate the optimal svalue of the fmincon into physical points
        boundary_point_of_interest(:,point_index) = projection_point;
        hyperplane_vector = (point_of_interest - projection_point)/...
                                norm(point_of_interest - projection_point);

        % The supporting hyperplane is given by: (P(x) - x_0)^T (x-x_0)<0
        Amatrix(point_index,:) = hyperplane_vector';
        bmatrix(point_index) = hyperplane_vector'*projection_point;
    end
    % Compute the polytopic overapproximation
    polytopic_overapproximation = Polyhedron('H', [Amatrix bmatrix]);
    % Compute the polytopic underapproximation
    polytopic_underapproximation = Polyhedron('V', boundary_point_of_interest');
end

function [c,ceq] = PrOccupySetConstraint(x, obstacle_mu, obstacle_sigma, ...
        obstacle_lower_bound, obstacle_upper_bound, log_alpha_value)
   % By proposition 3 and Lemma 5.1 --- Symmetric obstacle
   lb = (obstacle_lower_bound+x)';
   ub = (obstacle_upper_bound+x)';
   cdf_value = mvncdf(lb,ub, obstacle_mu', obstacle_sigma);
   c = log_alpha_value - log(max(cdf_value, eps));
   ceq = 0;
end

%%% Debugging helper: These lines when appropriately place help draw supporting
%%% hyperplanes at the boundary point
%% Setup
%step_size = 10;
%faraway_points_forward = zeros(2,no_of_points_on_the_circle);
%faraway_points_backward = zeros(2,no_of_points_on_the_circle);
%% At the end of the for-loop
%vector_on_the_line = cross([hyperplane_vector' 0], [0 0 1])/norm(cross([hyperplane_vector' 0], [0 0 1])); 
%faraway_points_forward(:,point_index) = projection_point+vector_on_the_line(1:2)'*step_size; 
%faraway_points_backward(:,point_index) = projection_point-vector_on_the_line(1:2)'*step_size;

%%% Debugging helper: Display the progress of fmincon (Add exitflag)
% constraint_value = nonlcon(projection_point);
% fprintf('Constraint: %1.3e | fmincon status: %d \n', constraint_value, exitflag);

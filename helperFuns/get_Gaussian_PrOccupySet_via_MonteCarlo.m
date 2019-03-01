%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name of the programmer: Abraham %
% Date: 2017-08-10                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Purpose
% Creates a grid and then computes the occupancy function over it via
% Monte-Carlo simulations
% Useful for creating a Monte-Carlo simulation-based superlevel set of the
% occupancy function

%% Inputs
% no_of_MC_points   : Number of Monte Carlo simulation particles
% no_of_grid_points : Number of grid points in the histogram (no. of bins)
% obstacle_mu       : Mean of the Gaussian FSRPD for the obstacle
% obstacle_sigma    : Covariance matrix of the Gaussian FSRPD for the obstacle
% obstacle_size     : Obstacle box size

%% Output
% xvec      : Grid points along x-axis
% yvec      : Grid points along y-axis
% frequency : Occupancy function estimate over the grid. For plotting purposes,
%               it is in the form frequency(y,x).

function [xvec,...
          yvec,...
          frequency,...
          other_values]=get_Gaussian_PrOccupySet_via_MonteCarlo(...
                                               no_of_MC_points,...
                                               no_of_grid_points,...
                                               obstacle_mu,...
                                               obstacle_sigma,...
                                               obstacle_size)

    boundary_distance_from_mean = 3*max(sqrt(eig(obstacle_sigma))) ...
                                    + 2 * obstacle_size;
    no_of_grid_points_x = no_of_grid_points;
    no_of_grid_points_y = no_of_grid_points;
    grid_min_x = obstacle_mu(1) - boundary_distance_from_mean;
    grid_max_x = obstacle_mu(1) + boundary_distance_from_mean;
    grid_min_y = obstacle_mu(2) - boundary_distance_from_mean;
    grid_max_y = obstacle_mu(2) + boundary_distance_from_mean;
    
    % Xvec and yvec creation
	xvec = linspace(grid_min_x,grid_max_x,no_of_grid_points_x);
    yvec = linspace(grid_min_y,grid_max_y,no_of_grid_points_y);

    % For index finding: Computation of step size
    xvec_step_size = xvec(2) - xvec(1);
    yvec_step_size = yvec(2) - yvec(1);

    % Translating the obstacle size into grid steps
    x_obstacle_size_in_grid_steps = round(obstacle_size/xvec_step_size);
    y_obstacle_size_in_grid_steps = round(obstacle_size/yvec_step_size);

    % Initialize the histogram matrix
    frequency = zeros(no_of_grid_points_y,no_of_grid_points_x);

    % Create mvn random samples
    mvn_samples = mvnrnd(obstacle_mu', obstacle_sigma, no_of_MC_points);
    fprintf('Created %d samples\n', no_of_MC_points);
    % Iterate over the samples
    for point_indx = 1: no_of_MC_points
        point_of_interest = mvn_samples(point_indx,:)';

        % Compute the index point and round it off to stay within xvec and yvec
        x_indx_point_orig = round((point_of_interest(1) - grid_min_x)/xvec_step_size)+1;
        y_indx_point_orig = round((point_of_interest(2) - grid_min_y)/yvec_step_size)+1;
        x_indx_point = min(max(x_indx_point_orig,1),no_of_grid_points_x);
        y_indx_point = min(max(y_indx_point_orig,1),no_of_grid_points_y);
        % Compute around the point_of_interest the box which is occupied the
        % obstacle
        x_iterates = max(x_indx_point - x_obstacle_size_in_grid_steps,1):min(x_indx_point + x_obstacle_size_in_grid_steps,no_of_grid_points_x);
        y_iterates = max(y_indx_point - y_obstacle_size_in_grid_steps,1):min(y_indx_point + y_obstacle_size_in_grid_steps,no_of_grid_points_y);

        % Increment the frequency in these states by 1
        frequency(y_iterates, x_iterates) = frequency(y_iterates, x_iterates) + 1;
    end
    disp('Completed the histogram');

    % Normalizing the histogram
    frequency = frequency/no_of_MC_points;

    other_values = [boundary_distance_from_mean,...
                    xvec_step_size,...
                    yvec_step_size,...
                    x_obstacle_size_in_grid_steps,...
                    y_obstacle_size_in_grid_steps];
end

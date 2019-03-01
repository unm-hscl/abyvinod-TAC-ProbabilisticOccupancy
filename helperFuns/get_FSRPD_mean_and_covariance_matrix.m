%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name of the programmer: Abraham %
% Date: 2017-08-10                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Purpose
% Computation of FSRPD via HSCC 2017, Proposition 1.

%% Notes
% For dynamics with additive Gaussian noise

%% Inputs
% ctrb_matrix               : Controllability matrix as given by the dynamics
% state_transition_matrix   : State transition matrix as given by the dynamics
% mu_disturbance            : Mean of the disturbance
% variance_disturbance      : Variance of the disturbance
% time_steps_taken          : Time for which FSRPD parameters must be computed
% no_of_noise_inputs        : No. of noise inputs in the dynamics (from
%                               controllability matrix this helps in estimating
%                               the time_horizon as well as size of B matrix)
% obstacle_init_state       : Initial state of the obstacle

%% Outputs
% mu_obstacle       : Mean of the FSRPD
% sigma_obstacle    : Sigma of the FSRPD

%% Tested by
% 1. test_get_FSRPD_mean_and_covariance_matrix.m

%% MODIFIED: ONLY HANDLES TIME_HORIZON COMPUTATION!

function [mu_obstacle, sigma_obstacle] = get_FSRPD_mean_and_covariance_matrix(...
                                           ctrb_matrix,...
                                           state_transition_matrix,...
                                           mu_disturbance,...
                                           variance_disturbance,...
                                           time_steps_taken,...
                                           obstacle_init_state)
    % TODO: assert time_steps_taken and ctrb_matrix are of the correct dimension
    obstacle_no_disturbance_reachable = state_transition_matrix * obstacle_init_state; 
    mu_obstacle = obstacle_no_disturbance_reachable +  ctrb_matrix * kron(ones(time_steps_taken,1), mu_disturbance);
    sigma_obstacle = ctrb_matrix * kron(eye(time_steps_taken), variance_disturbance) * ctrb_matrix';
    %if abs(det(sigma_obstacle))<eps
        %warning('Gaussian is degenerate in one dimension!')
    %end
end

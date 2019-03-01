%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name of the programmer: Abraham %
% Date: 2017-08-09                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Purpose
% Construct the controllability and state transition matrices for piecewise
% linear (PWL) system model of the unicycle

%% Notes
% The unicycle model is decomposed as a PWL system with its mode as the
%       turning rate \in {-W_max:W_step:W_max}.
% Using theta_init (initial heading) and mode_sequence_omega (the known 
%       sequence of turning rates), it is possible to predict theta for 
%       every time step.
% The sampling_time is the sampling time.

%% Inputs
% mode_sequence_omega   : Turning rate at each instant
% theta_init            : Initial heading
% sampling_time         : Sampling time

%% Outputs
% ctrb_matrix             : A matrix with two rows and 
%                           length(mode_sequence_omega)+1 columns
% state_transition_matrix : I_2, Identity matrix since without disturbance
%                           (velocity) the obstacle is not going to go
%                           anywhere

%% Tested by:
% 1. test_get_ctrb_and_state_transition_matrices_unicycle.m
% 2. test_get_FSRPD_mean_and_covariance_matrix.m

function [ctrb_matrix, state_transition_matrix] = get_ctrb_and_state_transition_matrices_unicycle(...
                                                    mode_sequence_omega,...
                                                    theta_init,...
                                                    sampling_time)
    % time_steps_taken is length(mode_sequence_omega) 
    ctrb_matrix = [];
    for i=1:length(mode_sequence_omega)+1
        if i == 1
            theta_current = theta_init;
        else
            theta_current = theta_current + mode_sequence_omega(i-1) * sampling_time;
        end
        B_matrix = sampling_time * [ cos(theta_current);
                                     sin(theta_current)];
        ctrb_matrix = [ctrb_matrix, B_matrix];        
    end
    state_transition_matrix = eye(2);

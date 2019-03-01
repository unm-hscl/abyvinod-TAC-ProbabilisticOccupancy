%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name of the programmer: Abraham %
% Date: 2017-08-11                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Purpose
% Get the probability of a given condensed_and_translated_mode_sequence_omega

%% Inputs
% condensed_and_translated_mode_sequence_omega  : 
% switching_probability_matrix                  : Markov transition probability
%                                                   governing the switching
% translated_initial_mode_omega                 : Translation of initial mode omega

%% Outputs
% probability_of_sequence   : Probability of occurence for the sequence

%% Notes
% Not vectorizable

function probability_of_sequence = get_probability_for_mode_sequence(...
                                        switching_sequence,...
                                        switching_probability_matrix,...
                                        indx_for_omega_init)
    if indx_for_omega_init ~= switching_sequence(1)
        % Initial mode does not match the given sequence => Improbable sequence
        probability_of_sequence = 0;
    else
        % Initial mode matches the given sequence => Probable sequence
        probability_of_sequence = 1;
        for time_for_switch= 1:length(switching_sequence)-1
            current_state = switching_sequence(time_for_switch);
            next_state = switching_sequence(time_for_switch + 1);
            % Extract the probability as the product of the correspondings
            % transition probabilities
            probability_of_sequence = probability_of_sequence *...
                                        switching_probability_matrix(...
                                                current_state,next_state);
        end        
    end

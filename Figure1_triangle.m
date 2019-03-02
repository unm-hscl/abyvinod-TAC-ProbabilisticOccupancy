% Computation previously took 21 seconds
% Obtain the triangle distribution-based FSR

clear;clc;close all;

% No plotting
cf_options.isPlot = false;

% Distribution of chance constraints --- make sure math is correct
joint_prob_thresh = 0.05;                   % \Delta
n_curr_dir = 20;                            % no_of_dir => Actual polytope can
                                            % have up to 2*n_curr_dir (up to
                                            % because CharFunTool may fail in
                                            % some directions 
% delta_i = \Delta/N/2   (/2 because we repeat both directions)
indiv_prob_thresh = joint_prob_thresh/n_curr_dir/2;

%% Problem parameters
sampling_time = 0.1;                        % Sampling time
heading_init = pi/10;                       % Initial heading 
init_location = zeros(2,1);                  
turning_rate_seq = [0,0,10,10,10,0,0,-5,-5,-5];
% Input space definition
vmax = 1;
avg_vel = 10;

%% Auxillary problem parameters --- derived from above
time_horizon = length(turning_rate_seq);
theta_seq = heading_init + cumsum(sampling_time * turning_rate_seq);
v_rv_pdf_obj = makedist('Triangular','a',-vmax,'b',0,'c',vmax);
v_rv = RandomVector('UserDefined', @(N) v_rv_pdf_obj.random(1,N));
input_vec = avg_vel * ones(time_horizon,1);

%% Dynamics
sys = LtvSystem('StateMatrix', @(t) eye(2), ...
                'InputMatrix', @(t) sampling_time * ...
                    [cos(theta_seq(t+1)); sin(theta_seq(t+1))], ...
                'InputSpace', Polyhedron('lb',-Inf,'ub',Inf), ...
                'DisturbanceMatrix', @(t) sampling_time * ...
                    [cos(theta_seq(t+1)); sin(theta_seq(t+1))], ...
                'Disturbance', v_rv);
            
%% Auxillary problem parameters continued
% Concatenated matrix            
[Z,~,G] = sys.getConcatMats(time_horizon);
% Direction vectors (Note that -pi/2 to pi/2 since we will check the symmetric
% direction during the computation
theta_vec = linspace(-pi/2,pi/2,n_curr_dir + 1)';
theta_vec = theta_vec(1:end-1);
all_curr_dir = [cos(theta_vec),sin(theta_vec)];

%% Plot setup
figure(1);
clf
hold on;
plot(Polyhedron('lb',-0.0001*ones(2,1),'ub',0.0001*ones(2,1)),'color','r')
plot(Polyhedron('lb',-0.0001*ones(2,1),'ub',0.0001*ones(2,1)),'color','g')
leg=legend('FSR set', '0.95-confidence set');
set(leg,'AutoUpdate','off','Location','SouthWest');
set(gca,'FontSize',25);
xlabel('x');
ylabel('y');
box on;

%% Monte carlo simulation
n_monte_carlo_sims = 1e5;
concat_state_realization = generateMonteCarloSims(n_monte_carlo_sims, sys, ...
    init_location, time_horizon, input_vec);
% Add to the plot the points from the Monte-Carlo simulation
for t=1:time_horizon+1
    scatter(concat_state_realization(2*t-1,:), ...
        concat_state_realization(2*t,:),30,'ko','filled')
end
axis equal;
drawnow;

%% Confidence region computation
timerVal = tic;
for t=1:time_horizon
    % Dynamics unrolled up to the point of interest
    ctrb_mat = G( 2*t-1: 2*t,1:t);
    control_ctrb_mat = ctrb_mat;
    state_mat = Z( 2*t-1: 2*t,:);

    % Unperturbed state
    x_unpert = state_mat * init_location + avg_vel.* sum(control_ctrb_mat,2);

    % Compute the FSR set
    relv_input_space = Polyhedron('lb',-vmax*ones(t,1),'ub',vmax*ones(t,1));
    fsr_set{t} = x_unpert + ctrb_mat * relv_input_space;

    % Compute the confidence region
    if rank(ctrb_mat) == 1      
        % Use the first row to remain --- using a heuristic that controllability
        % matrix first row serves as the direction along which the FSR set will
        % lie
        curr_dir = ctrb_mat(:,1)';

        % Get the CDF via CharFunTool
        zeroinput_zerostate_cf_struct = cf2DistGP( ...
            @(t) cf_TriangularSymmetric(t,1,1), [], ...
            [indiv_prob_thresh, 1 - indiv_prob_thresh],cf_options);            
        % MPT fails in enforcing equality constraints properly. So use vertex
        xlimits = x_unpert + curr_dir'.* zeroinput_zerostate_cf_struct.qf;
        plot(fsr_set{t}.V(:,1),fsr_set{t}.V(:,2),'r','LineWidth',5);
        plot(xlimits(1,:),xlimits(2,:),'g','LineWidth',5);
    else
        % Polytope defining matrices
        A = [];
        b = [];            
        for curr_dir_indx = 1:n_curr_dir    
            % Get the direction to be explored
            curr_dir = all_curr_dir(curr_dir_indx, :);

            % Get the 1xt vector that will multiply with the concatenated random
            % vector describing the velocities
            curr_dir_proj = curr_dir * ctrb_mat;
            
            % Get the CDF via CharFunTool
            zeroinput_zerostate_cf_struct = cf2DistGP( ...
                @(t) cf_TriangularSymmetric(t, curr_dir_proj), [], ...
                [indiv_prob_thresh, 1 - indiv_prob_thresh],cf_options);            
            % Define hyperplanes associated with this only if it returns reals
            if isreal(zeroinput_zerostate_cf_struct.qf)
                low_qf = min(zeroinput_zerostate_cf_struct.qf);
                high_qf = max(zeroinput_zerostate_cf_struct.qf);
                A = [A;
                     -curr_dir;
                     curr_dir];
                b = [b;
                     -curr_dir * x_unpert - low_qf;
                      curr_dir * x_unpert + high_qf];
            else
                % disp('Skipped');
            end
        end
        % Construct the polytope --- it is an overapproximation, so intersect it
        % with FSR set
        conf_region{t} = Polyhedron('A',A, 'b', b).intersect(fsr_set{t});
    end
end
toc(timerVal)

%% Validation via MC
for t= 3:time_horizon
    flag = conf_region{t}.contains(concat_state_realization(2*t+1:2*t+2,:));
    fprintf('Confidence level at t=%d: %1.3f\n', t, nnz(flag)/n_monte_carlo_sims);
    % Plot the FSR set
    plot(fsr_set{t},'alpha',0.5);
    % Plot the confidence region
    plot(conf_region{t},'alpha',0.4,'color','g');
end

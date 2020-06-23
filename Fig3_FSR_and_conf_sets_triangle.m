% Computation previously took 21 seconds
% Obtain the triangle distribution-based FSR

% We perform uniform risk allocation along n_curr_dir, so that collectively the 
% risk of the state NOT lying in the characterized polytope is below 0.05. We 
% use uniform risk allocation to split this joint chance constraint into 
% individual chance constraints using Boole's inequality. Consequently, along 
% each vector direction, we obtain two vertices of the polytopes via the 
% quantile function evaluated at indiv_prob_thresh and 1 - indiv_prob_thresh.
% Here, indiv_prob_thresh = (0.05/N * 2) with N direction vectors by
% uniform risk allocation (allocate risk to 2N hyperplanes). The resulting
% polytope (in V-Rep) is a guaranteed subset of the 0.95-confidence region.

clear;clc;close all;
diary off;
diary('figs/Figure3_console.txt');
diary on;
fprintf('\nNew experiment\n\n%s\n', datestr(now, 'DD mmmm, YYYY HH:MM:SS'));

% No plotting
cf_options.isPlot = false;

% Distribution of chance constraints --- make sure math is correct
n_monte_carlo_sims = 1e5;
skip_monte_carlo_plot = 1;                  % Skip plotting in scatter plot
joint_prob_thresh = 0.05;                   % \Delta
max_half_n_curr_dir = 20;                   % no_of_dir => Actual polytope can
                                            % have up to 2*n_curr_dir (up to
                                            % because CharFunTool may fail in
                                            % some directions 
% delta_i = \Delta/N/2   (/2 because we repeat both directions)
indiv_prob_thresh = joint_prob_thresh/max_half_n_curr_dir/2;

%% Problem parameters
sampling_time = 0.1;                           % Sampling time
initial_heading = pi/10;                       % Initial heading 
init_location = zeros(2,1);                  
turning_rate_seq = [0,0,10,10,10,0,0,-5,-5,-5];
% Input space definition
v_delta = 1;  % Should be a positive number
avg_vel = 10;

%% Auxillary problem parameters --- derived from above
time_horizon = length(turning_rate_seq);
v_rv_pdf_obj = makedist('Triangular','a',avg_vel-v_delta,'b', avg_vel, ...
    'c',avg_vel + v_delta);
v_rv = RandomVector('UserDefined', @(N) v_rv_pdf_obj.random(1,N));
% Direction vectors (Note that -pi/2 to pi/2 since we will check the
% opposite direction as well
theta_vec = linspace(-pi/2,pi/2,max_half_n_curr_dir + 1)';
theta_vec = theta_vec(1:end-1);
all_curr_dir = [cos(theta_vec),sin(theta_vec)];

%% Dynamics
% theta_seq = initial_heading + cumsum(sampling_time * [0,turning_rate_seq]);
% sys = LtvSystem('StateMatrix', @(t) eye(2), ...
%                 'DisturbanceMatrix', @(t) sampling_time * ...
%                     [cos(theta_seq(t+1)); sin(theta_seq(t+1))], ...
%                 'Disturbance', v_rv);
sys = getDubinsCarLtv('vel-dist', turning_rate_seq', initial_heading, ...
    sampling_time, v_rv);
% Concatenated matrix            
[Z,~,G] = sys.getConcatMats(time_horizon);

%% Plot setup
figure(1);
clf
hold on;
plot(Polyhedron('lb',-0.0001*ones(2,1),'ub',0.0001*ones(2,1)),'color','r')
plot(Polyhedron('lb',-0.0001*ones(2,1),'ub',0.0001*ones(2,1)),'color','g')
leg=legend('FSR set', '0.95-confidence set');
set(leg,'AutoUpdate','off','Location','SouthWest');
set(gca,'FontSize',40);
xlabel('x');
ylabel('y');
box on;

%% Confidence region computation
timerVal = tic;
for t=1:time_horizon
    % Dynamics unrolled up to the point of interest
    ctrb_mat = G( 2*t-1: 2*t,1:t);
    control_ctrb_mat = ctrb_mat;
    state_mat = Z( 2*t-1: 2*t,:);

    %% Compute the FSR set
    % Unperturbed state
    x_unpert = state_mat * init_location + avg_vel.* sum(control_ctrb_mat,2);
    % Zero mean version
    relv_dist_space = Polyhedron('lb',-v_delta*ones(t,1),'ub',v_delta*ones(t,1)); 
    % FSR set has t going from 0 to time_horizon (hence t+1)
    fsr_set{t+1} = x_unpert + ctrb_mat * relv_dist_space;

    % Compute the confidence region
    if rank(ctrb_mat) == 1      
        % Use the first row to remain --- the first column of the 
        % controllability matrix is the direction along which the FSR set will 
        % lie. Other columns must be scaled versions of this.
        curr_dir = ctrb_mat(:,1)';

        % Get the confidence interval via CharFunTool by inverting the CDF
        % at confidence threshold beta/2.
        zeroinput_zerostate_cf_struct = cf2DistGP( ...
            @(t) cf_TriangularSymmetric(t,1,1), [], ...
            [indiv_prob_thresh, 1 - indiv_prob_thresh],cf_options);            
        
        % MPT fails in enforcing equality constraints properly. So compute
        % the vertices directly
        xlimits = x_unpert + curr_dir'.* zeroinput_zerostate_cf_struct.qf;
        % conf_region has t going from 0 to time_horizon (hence t+1)
        conf_region{t+1} = Polyhedron('V', xlimits');
        
        % Plot them as lines, since polytope version is just a thin line
        plot(fsr_set{t+1}.V(:,1), fsr_set{t+1}.V(:,2),'r','LineWidth', 5);
        plot(xlimits(1,:),xlimits(2,:),'g','LineWidth', 5);
        drawnow;
    else
        % Polytope defining matrices
        A = [];
        b = [];            
        for curr_dir_indx = 1:max_half_n_curr_dir    
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
%                 disp(zeroinput_zerostate_cf_struct.qf)
                low_qf = min(zeroinput_zerostate_cf_struct.qf);
                high_qf = max(zeroinput_zerostate_cf_struct.qf);
                A = [A;
                     -curr_dir;
                     curr_dir];
                b = [b;
                     -curr_dir * x_unpert - low_qf;
                      curr_dir * x_unpert + high_qf];
            else
%                 fprintf(['Skipped a direction (t=%d, dir indx = %d) since ', ...
%                     'imaginary limits found!\n'], t, curr_dir_indx);
            end
        end
        % Construct the polytope --- it is an overapproximation, so intersect it
        % with FSR set
        conf_region{t+1} = Polyhedron('A',A, 'b', b).intersect(fsr_set{t+1});
    end
end
toc(timerVal)


%% Monte carlo simulation
concat_state_realization = generateMonteCarloSims(n_monte_carlo_sims, sys, ...
    init_location, time_horizon);
% Add to the plot the points from the Monte-Carlo simulation
for t=1:time_horizon+1
    scatter(concat_state_realization(2*t-1,1:skip_monte_carlo_plot:end), ...
            concat_state_realization(2*t,1:skip_monte_carlo_plot:end), ...
            10, 'ko', 'filled')
end
axis tight;
drawnow;

%% Validation via MC for sets at t\in[4,T], T is the time horizon
% t={1,2,3} are one-dimensional and plotted already
for t= 4:time_horizon
    % concat_state_realizations has (T+1) rows to include the initial state
    mc_realizations = concat_state_realization(2*t+1:2*t+2,:);
    % conf_regions have T+1 polytopes
    flag = conf_region{t+1}.contains(mc_realizations);
    % fsr_set have T+1 polytopes
    plot(fsr_set{t+1},'alpha',0.5);
    % Plot the confidence region
    plot(conf_region{t+1},'alpha',0.4,'color','g');
    fprintf('Confidence level at t=%d: %1.3f\n', t, ...
        nnz(flag)/n_monte_carlo_sims);    
end
axis equal;
savefig(gcf, 'figs/Unicycle_FSR_and_95_conf_region.fig', 'compact');
saveas(gcf, 'figs/Unicycle_FSR_and_95_conf_region.png');
diary off;
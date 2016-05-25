
function  clf = scaling_chain()
                                
%         '''
%         m: integrator chain order
%         n: dimension of output space = number of control inputs
%         A(m*n,m*n): State matrix for the given integrator chain
%         B(m*n,n): Control matrix for the given integrator chain
%         obstacles: list of corners of the obstacles
%         goal: list of corners of the goal location
%         X0: initial state vector
%         Realizability: Set to False if there problem may not be solved for some reason, terminates the problem
%         '''
prob = problem_init;
prob = prob.init_prob();
repeat = 'y';
while repeat == 'y'
clf = CLF_CBF_QP;
clf.n = prob.output_dim;
clf.m = prob.number_integrators;
clf.goal_num = prob.goal_num;
clf.obs_num = prob.obstacle_num;
clf = clf.initialize();

targets = zeros(clf.n,prob.goal_num);
obs_size = zeros(clf.n*2,prob.obstacle_num);
% '''
% Get extremums of each obstacle and store as a list in "box_size"
% Define obstacles
% '''
for i = 1:prob.obstacle_num
    obs = prob.obstacles(:,:,i);
    obs_min = min(obs);
    obs_max = max(obs);
    obs_temp = [];
    for j = 1:clf.n
        obs_temp = [obs_temp,obs_min(j)];
        obs_temp = [obs_temp,obs_max(j)];
    end
    obs_size(:,i) = obs_temp;
end
% '''
% Get feasible target locations and store as a list in "targets"
% Even if one location is not feasible the the problem in not "Realizable"
% '''
clf.box_size = obs_size;
goal_size = zeros(clf.n*2,prob.goal_num);

for i = 1:prob.goal_num
    target_polytopeV = prob.goals(:,:,i);
    target = random_goal(target_polytopeV,clf.box_size);
    goal_min = min(target_polytopeV);
    goal_max = max(target_polytopeV);
    goal_temp = [];
    for j = 1:clf.n
        goal_temp = [goal_temp,goal_min(j)];
        goal_temp = [goal_temp,goal_max(j)];
    end
    goal_size(:,i) = goal_temp;
    if isempty(target)
        fprintf('initial position/goal in collision with box.. random sampling failed\n')
        clf.Realizability = false;
        break;
    else
        targets(:,i) = target;
    end
end
clf.goal_size = goal_size;
clf.X0 = prob.X0;
clf.U_min = prob.U_min;
clf.U_max = prob.U_max;
% '''
% Construct the CBF vector for outter bounds Y
% - Get the bounds
% - Find the ellipse that inscribes the bounds which is the CBF: B(x)
% - Generate the CBF state vector (CBF_Y): [B(x); dB(x); .. d^(m-1)B(x)]
% - Generate the feedback gain matrix K_Y which ensures we remain in bounding box y<=0
% '''
A_cbf = diag(ones(clf.m-1,1),1);
B_cbf = [zeros(clf.m-1,1);1];
p = -1:-1:-clf.m;
Ks = [];
cbf_p(clf.obs_num) = cbf_prms;
K_Ys = [];

% %get Y limits
clf.Y_min = prob.Y_min;
clf.Y_max = prob.Y_max;
if clf.Realizability
    for i =1:clf.n
        eta0_Y = clf.X0(i:clf.n:end);
        eta0_Y(1) = eta0_Y(1)-clf.Y_min(i);
%         fprintf('B_Y0:', eta0_Y);
        if eta0_Y(1) < 0
            fprintf('start location outside bounding box\n');
%             fprintf('initial state:', clf.X0);
            clf.Realizability = false;
        end
        if clf.m == 1
            K_Y = 1;
        else
            K_Y = CBF_pole_design(eta0_Y,A_cbf,B_cbf,p);
%             fprintf('generated K_Y:',K_Y)
            
        end
        K_Ys = [K_Ys;K_Y];
        eta0_Y = clf.X0(i:clf.n:end);
        eta0_Y(1) = eta0_Y(1)-clf.Y_max(i);
%         fprintf('B_Y0:', eta0_Y);
        if eta0_Y(1) > 0
            fprintf('start location outside bounding box\n');
%             fprintf('initial state:', clf.X0);
            clf.Realizability = false;
        end
        if clf.m == 1
            K_Y = 1;
        else
            K_Y = CBF_pole_design_Y(eta0_Y,A_cbf,B_cbf,p);
%             fprintf('generated K_Y:',K_Y)
        end
        K_Ys = [K_Ys;K_Y];
        
    end
    clf.K_Ys = K_Ys;
    % '''
    % Construct the CBF vector for each obstacle
    % - Get the bounds
    % - Find the ellipse that circumscribes the bounds which is the CBF: B(x)
    % - Generate the CBF state vector (cbf): [B(x); dB(x); .. d^(m-1)B(x)]
    % - Generate the feedback gain matrix K which ensures we remain in bounding box y<=0
    % - Store each K and cbf in a list of "Ks" and "cbf_prms"
    % '''
    for i = 1:clf.obs_num
        box = obs_size(:,i);
        ellipse = ellips_obs_constr(box,targets,clf.n,0,clf.X0);
        if ~isempty(ellipse.center)
            cbf = cbf_p(i).cbf_gen(clf.m,clf.n,ellipse);
            cbf_p(i) = cbf;
            cbf_sub = cbf_p(i).subs_cbf(clf.X0);
            eta0 = cbf_sub.eta;
            fprintf('B0: %d\n',eta0(1,1));
            if eta0(1,1) <= 0
                if chk_collision(clf.X0(1:clf.n),box)
                    fprintf('start location in obstacle\n');
                    clf.Realizability = false;
                end
            end
            if clf.m == 1
                K = 1;
            else
                K = CBF_pole_design(eta0,A_cbf,B_cbf,p);
                fprintf('\ngenerated K: ');
                disp(K);
            end
            Ks = [Ks;K];
            
        else
            trail_ended = true;
            clf.Realizability = false;
            fprintf('Trail ended because goal is very close to obstacle\n');
        end
    end
    clf.Ks = Ks;
    clf.cbf_p = cbf_p;
end
% '''
% Initialize the control and check till the problem converges (error<0.01) or till the time ellapses
% If problem converges try to reach the next goal region
% '''
goal_num = prob.goal_num;
current = 1;
t0 = 0;
tf = 100;
traj = [];
traj_t = [];
if clf.Realizability
    while (goal_num-current) ~= -1
        clf.goal = targets(:,current);
        clf.target_state = [clf.goal; zeros((clf.m-1)*clf.n,1)];
        clf.vel = (clf.goal - clf.X0(1:clf.n))./clf.T;  
        options = odeset('Events',@event_func);
        [t,x_dynam] = ode45(@f_nonlin,[t0 tf],clf.X0,options,clf);
        clf.X0 = x_dynam(end,:);
        t0 = t(end);
        fprintf('Reached Goal %d\n', current+1);
        current = current+ 1;
        traj = [traj; x_dynam];
        traj_t = [traj_t;t];
    end
else
    fprintf('trail ended because implementation is not realizable\n');
end
repeat = input('Do you want to repeat the problem?(y/n)','s');
end
end

function target = random_goal(goal,box_size)
% '''
% goal: goal region with corner points of polytope
% box_size: contains extreme values of the polytope per dimension
% sample random points in the goal region and find a "target" which is not in collision
% '''
goal_min = min(goal);
goal_max = max(goal);
target = [];
n = length(goal_min);
ran_num = 20^n;
rand_goals = repmat(goal_min,ran_num,1) + repmat((goal_max-goal_min),ran_num,1).*rand(ran_num,n);
for j = 1:ran_num
    g = rand_goals(j,:);
%     for i = 1:n
%         g(i) = goal_min(i) + g(i)*(goal_max(i)-goal_min(i));
%     end
    if ~chk_collision(g,box_size)
        target = g;
        return;
    end
end
                                                                 
end


                                                        
function K = CBF_pole_design(eta0,A,B,p)
% '''
% eta0 is the initial value of B,bB,....d^(m-1)B got from X0
%     A and B are the dynamics state and control matrices
% p contains m pole locations which are mostly chosen to be negative
% output K is the feedback chosen such that the poles ensure that the state eta will always remain outside(positive) the polytope
% '''
y0 = CBF_output(eta0,p);
% % Modify poles until all elements of y0 >= 0
while ~all(y0>=0)
        p = 1.1*p;  %% then modify pole location
        y0 = CBF_output(eta0,p);
end
K = place(A, B, p);
end       

function K = CBF_pole_design_Y(eta0,A,B,p)
% '''
% eta0 is the initial value of B,bB,....d^(m-1)B got from X0
%     A and B are the dynamics state and control matrices
% p contains m pole locations which are mostly chosen to be negative
% output K is the feedback chosen such that the poles ensure that the state eta will always remain inside(negative) the polytope
% '''
y0 = CBF_output(eta0,p);
% % Modify poles until all elements of y0 >= 0
while ~all(y0<=0)
        p = 1.1*p;  %% then modify pole location
        y0 = CBF_output(eta0,p);
end
K = place(A, B, p);
end
                                
function y = CBF_output(eta,p)
% '''
% eta is a vector of dimension m containing values of B,dB,...,d^(m-1)B
%     p contains m pole locations which are mostly chosen to be negative
% output 	y(0) = eta(0) = B
% y(r) = (D-p1)(D-p2)...(D-pr)B
% '''
m = length(eta);
y = zeros(m-1,1);
y(1) = eta(1);
for i = 1:m-1
    coeff = poly(p(1:i));
    y(i+1) = coeff(end:-1:1)*eta(1:i+1);
end
end
        
function ellipse = ellips_obs_constr(box_size, goal, n, obs_flag,X0)
% '''
% Generate an ellipsoid-shaped obstacle to approximate a box
% box_size: a list containing the minimum and maximum of each dimension
% goal:     a list containing each goal position
% n:    the dimension of the position space
% obs_flag:	the flag determine whether the box should be inscribed(1) or circumscribed(0)
% '''
ellipse = hyper_ellipse;
goal_num = size(goal,2);
goal(:,goal_num+1) = X0(1:n);
goal_num = goal_num+1;
if obs_flag == 0
    p = 1.0;
else
    p = 1.0;               ...  % the initial degree of this super ellipsoid
end
max_p = 2^5;              ...  % the maximal number of p's degree
box_center = (box_size(1:2:end)+box_size(2:2:end))./2;
dist = box_center - box_size(1:2:end);
flag = false;
while flag ~= true
    p = p * 2; ...   % double the power
        if p > max_p
            % if the power raises too high, the ellipsoid approximation is not good enough
            % stop it to prevent infinite loop
            fprintf('Cannot approximate this box with degree less than %d\n',max_p);
            return
        end
        a = dist.^p;
        flag = true;
        
        for j = 1:goal_num
            exo = box_center - goal(j);
            g = (exo.^p)./a;
            g = sum(g);
            if g <= n && obs_flag == 0
                % one goal lies within the obstacle
                flag = false;
                break;
            end
            if g > 1 && obs_flag == 1
                %% one goal lies outside the obstacle
                flag = false;
                break;
            end
        end
end

% scale the vector a for better numerical stability
% and store the values to geom
scale = min(a);
a = a./scale;
if obs_flag == 0
    b = n * scale;
else
    b = scale;
end

ellipse.center = box_center;
ellipse.axis = a;
ellipse.const = b;
ellipse.degree = p;
end



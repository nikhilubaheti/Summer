classdef problem_init
    properties (Access = public)
        output_dim_range = [2 2];
        number_integrators_range = [1 4];
        goal_range = [1 10];
        obstacle_range = [1 4];
        t0 = 0;
        tf = 3;
        output_dim = 1;
        number_integrators = 1;
        U_min = []; 
        U_max = [];
        Y_min = [];
        Y_max = [];
        U_box = [];
        Y_box = [];
        goal_num = 1;
        obstacle_num = 1;
        goals = [];
        obstacles = [];
        X0 = [];

    end
    methods
        function prob = init_prob(prob)
            prob.output_dim = randi(prob.output_dim_range);
            prob.number_integrators = randi(prob.number_integrators_range);
            prob.U_min = -ones(prob.output_dim,1);
            prob.U_max = ones(prob.output_dim,1);
            prob.Y_min = -10*ones(prob.output_dim,1);
            prob.Y_max = 10*ones(prob.output_dim,1);
            prob.U_box = generate_box(prob.U_min,prob.U_max)';
            prob.Y_box = generate_box(prob.Y_min,prob.Y_max)';
            prob.goal_num = randi(prob.goal_range);
            prob.obstacle_num = randi(prob.obstacle_range);
            prob.goals = generate_boxes(prob.goal_num,prob.Y_min,prob.Y_max);
            prob.obstacles = generate_boxes(prob.obstacle_num,prob.Y_min,prob.Y_max);
            prob.X0 = [prob.Y_min+(prob.Y_max-prob.Y_min).*rand(prob.output_dim,1); zeros((prob.output_dim)*(prob.number_integrators-1),1)];
        end
        
        function obs_size = obstacle_extremums(prob)
            % '''
            % Get extremums of each obstacle and store in "obs_size"
            % '''
            obs_size = zeros(prob.output_dim*2,prob.obstacle_num);
            for i = 1:prob.obstacle_num
                obs = prob.obstacles(:,:,i);
                obs_min = min(obs);
                obs_max = max(obs);
                for j = 1:prob.output_dim
                    obs_size(2*j-1,i) = obs_min(j);
                    obs_size(2*j,i) = obs_max(j);
                end
            end
        end
        
        function [targets,goal_size,Realizability] = target_points(prob,clf,obs_size)
            % '''
            % Get feasible target locations and store in "targets"
            % Even if one location is not feasible the the problem in not "Realizable"
            % Get extremums of each goal and store in "goal_size"
            % '''
            
            targets = zeros(prob.output_dim,prob.goal_num);
            goal_size = zeros(prob.output_dim*2,prob.goal_num);
            Realizability = true;
            for i = 1:prob.goal_num
                target_polytopeV = prob.goals(:,:,i);
                
                goal_min = min(target_polytopeV);
                goal_max = max(target_polytopeV);
                for j = 1:prob.output_dim
                    goal_size(2*j-1,i) = goal_min(j);
                    goal_size(2*j,i) = goal_max(j);
                end
                target = goal_intersection_obs(goal_size(:,i),obs_size,clf);
%                 target = random_goal(target_polytopeV,obs_size);
                if isempty(target)
%                     fprintf('initial position/goal in collision with box.. random sampling failed\n');
                    fprintf('goal in collision with box\n');
                    Realizability = false;
                    break;
                else
                    targets(:,i) = target;
                end
            end
        end
        
        function [K_Ys,Realizability] = state_bounds(prob)
            % '''
            % Construct the CBF vector for outter bounds Y
            % - Get the bounds
            % - Generate the CBF state vector (eta0_Y): [B(x); dB(x); .. d^(m-1)B(x)]
            % - B(x) = Hx-K: K is the bounds H is identity
            % - Generate the feedback gain matrix K_Y which ensures we remain in bounding box y<=0
            % '''
            Realizability = true;
            A_cbf = diag(ones(prob.number_integrators-1,1),1);
            B_cbf = [zeros(prob.number_integrators-1,1);1];
            p = -1:-1:-prob.number_integrators;
            K_Ys = zeros(prob.output_dim*2,prob.number_integrators);
            for i =1:prob.output_dim
                eta0_Y = prob.X0(i:prob.output_dim:end);
                eta0_Y(1) = eta0_Y(1)-prob.Y_min(i);
                if eta0_Y(1) < 0
                    fprintf('start location outside bounding box\n');
                    Realizability = false;
                end
                if prob.number_integrators == 1
                    K_Y = 1;
                else
                    K_Y = CBF_pole_design(eta0_Y,A_cbf,B_cbf,p);
                    
                end
                K_Ys(2*i-1,:) = K_Y;
                eta0_Y = prob.X0(i:prob.output_dim:end);
                eta0_Y(1) = eta0_Y(1)-prob.Y_max(i);
                if eta0_Y(1) > 0
                    fprintf('start location outside bounding box\n');
                    Realizability = false;
                end
                if prob.number_integrators == 1
                    K_Y = 1;
                else
                    K_Y = CBF_pole_design_Y(eta0_Y,A_cbf,B_cbf,p);
                end
                K_Ys(2*i,:) = K_Y;
                
            end
        end
        
        function [Ks,cbf_p,Realizability,ellipses] = obstacle_bounds(prob,clf,obs_size,targets)
            % '''
            % Construct the CBF vector for each obstacle
            % - Find the ellipse that circumscribes the bounds which is the CBF: B(x)
            % - Generate the CBF state vector (cbf): [B(x); dB(x); .. d^(m-1)B(x)]
            % - Generate the feedback gain matrix K which ensures we remain in bounding box y<=0
            % - Store each K and cbf in "Ks" and "cbf_prms"
            % '''
            Ks = zeros(prob.obstacle_num,prob.number_integrators);
            cbf_p(prob.obstacle_num) = cbf_prms;
            ellipses(prob.obstacle_num) = hyper_ellipse;
            A_cbf = diag(ones(prob.number_integrators-1,1),1);
            B_cbf = [zeros(prob.number_integrators-1,1);1];
            p = -1:-1:-prob.number_integrators;
            Realizability = true;
            for i = 1:prob.obstacle_num
                box = obs_size(:,i);
                ellipse = hyper_ellipse;
                ellipse = ellipse.ellips_obs_constr(clf,box,targets,prob.output_dim,0,prob.X0);
                if ~isempty(ellipse.center)
                    ellipses(i) = ellipse;
                    cbf = cbf_p(i).cbf_gen(prob.number_integrators,prob.output_dim,ellipse);
                    cbf_p(i) = cbf;
                    cbf_sub = cbf_p(i).subs_cbf(prob.X0);
                    eta0 = cbf_sub.eta;
                    fprintf('B0: %d\n',eta0(1,1));
                    if eta0(1,1) <= 0
                        if chk_collision(prob.X0(1:prob.output_dim),box)
                            fprintf('start location in obstacle\n');
                            Realizability = false;
                        else
                            eta0(1,1) = abs(eta0(1,1)); %X0 is not in obs but in the ellipse region so replace -ve eta by positive
                        end
                    end
                    if prob.number_integrators == 1
                        K = 1;
                    else
                        K = CBF_pole_design(eta0,A_cbf,B_cbf,p);
                        fprintf('\ngenerated K: ');
                        disp(K);
                    end
                    Ks(i,:) = K;
                    
                else
                    Realizability = false;
                    fprintf('Trail ended because goal is very close to obstacle\n');
                end
            end
        end
    end
end


function boxes = generate_boxes(box_num,dim_min,dim_max)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input:
%     box_num: Generates box_num number of boxes where 
%     dim_min(n,1): contains the minimum values of the box planes for each dimension
%     dim_max(n,1): contains the maximum values of the box planes for each dimension
% Output: 
%     boxes(2^n,n,box_num): contains the corners of each box
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
output_dim = length(dim_min);
boxes = zeros(2^output_dim,output_dim,box_num);
for i = 1:box_num
    dim_rand = repmat(dim_min,1,2)+repmat((dim_max-dim_min),1,2).*rand(output_dim,2);
    goal_min = min(dim_rand,[],2);
    goal_max = max(dim_rand,[],2);
    boxes(:,:,i) = generate_box(goal_min,goal_max)';
end
end


function box = generate_box(dim_min,dim_max)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input:
%     dim_min(n,1): contains the minimum values of the box planes for each dimension
%     dim_max(n,1): contains the maximum values of the box planes for each dimension
% Output: 
%     box(n,2^n): contains the corners of the box
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = length(dim_min);
if n == 1
    box = [dim_min dim_max];
else
    box = combvec([dim_min(1) dim_max(1)],[dim_min(2) dim_max(2)]);
    for i = 3:n
        box = combvec(box,[dim_min(i) dim_max(i)]);
    end
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

function target = goal_intersection_obs(goal,obs_size,clf)
% '''
% Input: 
%   goal(2*n,1): goal region with extremum points of polytope
%   obs_size(2*n,obs_num): contains extreme values of the polytope per dimension
% Output:
%   target(n,1): A point in goal region not in present in the obstacle region.
%                     If entire goal region is in obstacle then the target is [].
% Tuning: 
%     choose_corners: true => target point is chosen as a corner point if the goal region is in obs
%                     but if all 4 corners are in obs then choose mean
%                     false => target point is chosen as the mean of the goal region which is not is obs
% '''
choose_corner = clf.commands.choosing_targets_corners;
target = [];
goal_sub_regions = goal;
n = length(goal)/2;
obs_length = size(obs_size,2);
region_no = 1; %keeps track of the current region that is being tested. Helps to check if the region tested has been
                   %check for all obstacles
obs_num = 1; %keeps track of current obstacle region being checked
obs_chked = [];%ensures that the current goal region has been checked upto this obstacle number
while obs_num <= obs_length
    current_goal = goal_sub_regions(:,region_no); 
    box = obs_size(:,obs_num);
    obs_chked(region_no) = obs_num; %stores the obs_num that the goal_region has been chked for
    obs_sub_region = []; %contains the sub goal regions for each obstacle
    Flag_dim = false(n,1); %check for every dimension if it lies in obstacle region
    for i = 1:n
        if box(2*i-1)>current_goal(2*i-1) && box(2*i-1)<current_goal(2*i) %check if ob_min lies within goal region; left region exists
            new_region = current_goal;
            new_region(2*i) = box(2*i-1); %the box must be generated on the left (min -> dim_obs)
            obs_sub_region = [obs_sub_region new_region];
            Flag_dim(i) = true;
        end
        if box(2*i)>current_goal(2*i-1) && box(2*i)<current_goal(2*i) %check if ob_max lies within goal region; right region exists
            new_region = current_goal;
            new_region(2*i-1) = box(2*i); %the box must be generated on the right (dim_obs -> max)
            obs_sub_region = [obs_sub_region new_region];
            Flag_dim(i) = true;
        end
        
        if Flag_dim(i) == false && (box(2*i-1)>current_goal(2*i-1) || box(2*i)<current_goal(2*i)) %check if obs_min>g_max or obs_max<g_min;
            obs_sub_region = []; %reset because the obs region is not in obstacle
            break;
        else
            Flag_dim(i) = true;
        end
    end
    if all(Flag_dim) %obstacle lies within goal region
        goal_sub_regions = [goal_sub_regions(:,1:region_no) obs_sub_region goal_sub_regions(:,region_no+1:end)];%modify goal regions with new goal_regions
        obs_chked = [obs_chked(1:region_no); obs_num*ones(size(obs_sub_region,2),1);obs_chked(region_no+1:end)];
        region_no = region_no+1; %move to the next goal region since current region is in collision
        
        if region_no > size(goal_sub_regions,2)
            warning('Goal region in obstacle');
%             target = [];
%             figure;
%             title(n);
%             hold on;
%             plot_boxes(n, [obs_size(1:2:end,:); obs_size(2:2:end,:)],'red',0.3,true) ;
%             plot_boxes(n, [goal_sub_regions(1:2:end,:); goal_sub_regions(2:2:end,:)],'green',0.3,true) ;
%             hold off;
%             pause;
%             close all;
            return;
        else 
            obs_num = obs_chked(region_no);
        end
    end
obs_num = obs_num+1;
end
current_goal = goal_sub_regions(:,region_no);  % update again to get latest region
if clf.commands.disp_goal_free
    figure;
    title(n);
    hold on;
    plot_boxes(n, [obs_size(1:2:end,:); obs_size(2:2:end,:)],'red',0.3,true) ;
    plot_boxes(n, [current_goal(1:2:end); current_goal(2:2:end)],'green',0.3,true) ;
    hold off;
end
if choose_corner && region_no > 1
    goal_box = generate_box(current_goal(1:2:end),current_goal(2:2:end));
    corner_flag = false;
    for i = 1:2^n
        if ~chk_collision(goal_box(:,i),obs_size)
            target = (goal_box(:,i)+mean(goal_box,2))'/2;
            corner_flag = true;
            break;
        end
    end
    if ~corner_flag
        target = (current_goal(1:2:end) + current_goal(2:2:end))/2';
    end
else
    target = (current_goal(1:2:end) + current_goal(2:2:end))/2';
    target = target';
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
while ~all(y0(2:end)>=0)
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
while ~all(y0(2:end)<=0) %%TO_DO verify if this if fine
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


    
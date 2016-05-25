classdef problem_init
    properties (Access = public)
        output_dim_range = [2 3];
        number_integrators_range = [1 4];
        goal_range = [1 10];
        obstacle_range = [1 4];
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
    end
end


function boxes = generate_boxes(box_num,Y_min,Y_max)
output_dim = length(Y_min);
boxes = zeros(2^output_dim,output_dim,box_num);
for i = 1:box_num
    dim_rand = repmat(Y_min,1,2)+repmat((Y_max-Y_min),1,2).*rand(output_dim,2);
    goal_min = min(dim_rand,[],2);
    goal_max = max(dim_rand,[],2);
    boxes(:,:,i) = generate_box(goal_min,goal_max)';
end
end


function box = generate_box(dim_min,dim_max)
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
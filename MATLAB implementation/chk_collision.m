function flag = chk_collision(goal,box_size)
% '''
% goal: goal position
% box_size: contains extreme values of the polytope per dimension
% check if goal location is in collision
% Check if for all dimensions the goal region lies in any obstacle(box)
% '''
flag = false;
n = length(goal);
for j = 1:size(box_size,2)
    Flag_dim = false(n,1);
    box = box_size(:,j);
        for i =1:n
            if goal(i)>=box(2*i-1) && goal(i)<=box(2*i)
                Flag_dim(i) = true;
            end
        end
        if all(Flag_dim)
            flag = true;
            return;
        end
end
end
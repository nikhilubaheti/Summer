function [value,isterminal,direction] = event_func(t,x,clf)
% clf.Reached_Goal = chk_collision(x(1:clf.n),clf.goal_size);
if clf.Reached_Goal
    value = 0;
else
    value = norm(clf.target_state - x);
end
isterminal = 1;
direction = 0;
end
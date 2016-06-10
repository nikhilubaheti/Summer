function display_trail(clf,t_reach_goal)
    fprintf('output_dim: %d\n',clf.n);
    fprintf('number of integrators:%d\n',clf.m);
    for i = 1:length(t_reach_goal)
        fprintf('reached goal %d ',i);
        fprintf('at %ds\n',t_reach_goal(i));
    end
end
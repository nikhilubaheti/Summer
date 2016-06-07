% Initialize the problem and the controller
% profile on -history
clc;
close all;
repeat = input('Do you want to repeat problem?(y/n)','s');
if repeat == 'y'
    load('prob.mat');
    prob.tf = 20;
%     prob.X0 = [-9.9 9.9 -9.9]';
else
    prob = problem_init;
    prob = prob.init_prob();
end
obs_size = prob.obstacle_extremums();
[targets,goal_size,Realizability] = prob.target_points(obs_size);
K_Ys = [];
Ks =[];
cbf_p = [];
Realizability = false;
if Realizability
    [K_Ys,Realizability] = prob.state_bounds();
end
if Realizability
    [Ks,cbf_p,Realizability] = prob.obstacle_bounds(obs_size,targets);
end
current = 1;
traj_x = [];
traj_u = [];
traj_t = [];
t0 = prob.t0;
t_last = prob.t0;
clf = CLF_CBF_QP;
clf = clf.initialize(prob,obs_size,goal_size,K_Ys,Ks,cbf_p,Realizability);

if clf.Realizability
    while (prob.goal_num-current) ~= -1 && t_last < prob.tf
        clf.goal = targets(:,current);
        clf.current_goal = current;
        clf.target_state = [clf.goal; zeros((clf.m-1)*clf.n,1)];
        clf.vel = (clf.goal - clf.X0(1:clf.n))./clf.T;
        options = odeset('Events',@event_func);
        [t,x_dynam] = ode15s(@f_nonlin,[t0 prob.tf-t_last],clf.X0,options,clf);
        clf.X0 = x_dynam(end,:)';
        dx = diff(x_dynam)./repmat(diff(t),1,clf.m*clf.n);
        u = dx(:,clf.n*(clf.m-1)+1:end);
        traj_u = [traj_u;u];
        traj_x = [traj_x; x_dynam(:,1:clf.n)];
        traj_t = [traj_t;t_last+t];
        t_last = traj_t(end);
        if t_last < prob.tf
            fprintf('Reached Goal %d\n', current);
            current = current+ 1;
        end
    end
    clf.X0 = prob.X0;
    plot_trials(clf,traj_u,traj_x);
else
    fprintf('trail ended because implementation is not realizable\n');
end
save('prob.mat','prob');
% p = profile('info');
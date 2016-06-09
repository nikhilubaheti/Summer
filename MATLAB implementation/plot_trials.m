function plot_trials(clf,traj_u,traj_y,traj_t,targets)
if(clf.n == 2)
    plot_trial(2,clf,traj_u,traj_y,traj_t,targets) ;
elseif(clf.n == 3)
    plot_trial(3,clf,traj_u,traj_y,traj_t,targets) ;
end
end


function plot_trial(n_dim,clf,u,y,t,targets)
figure ;
title('u plotted are approximation got from x_dynam');
for i = 1:clf.n
    subplot(clf.n,1,i);
    plot(t(1:size(u,1)),u(:,i)) ; grid on ; title(['u(' num2str(i) ')']);
end
    
    
figure ;
if(n_dim==2)
    plot(y(:,1), y(:,2), 'LineWidth', 5) ; hold on ;
    plot(clf.X0(1), clf.X0(2), 'k*') ;
    plot(targets(1,:), targets(2,:), 'b*') ;
else
    plot3(y(:,1), y(:,2), y(:,3), 'LineWidth', 5) ; hold on ;
    plot3(clf.X0(1), clf.X0(2), clf.X0(3), 'k*') ;
    plot3(targets(1,:), targets(2,:), targets(3,:), 'b*') ;
    zlabel('z') ;
end
xlabel('x'); ylabel('y') ; grid on ;
    
plot_boxes(n_dim, [clf.box_size(1:2:end,:); clf.box_size(2:2:end,:)],'red',0.3,false) ;
plot_boxes(n_dim, [clf.goal_size(1:2:end,1:clf.current_goal); clf.goal_size(2:2:end,1:clf.current_goal)],'green',0.3,true) ;
plot_boxes(n_dim, [clf.Y_min; clf.Y_max],'blue',0.01,false) ;
plot_ellipses(n_dim, clf.ellipses);
end

function plot_ellipses(n,ellipses)
num = length(ellipses);
if n == 2
    for i = 1:num
        if isempty(ellipses(i).center)
            break;
        end
        plot_se(ellipses(i).center,(ellipses(i).axis*ellipses(i).const).^(1/ellipses(i).degree),ellipses(i).degree);
        
    end
end
end
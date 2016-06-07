function plot_trials(clf,traj_u,traj_y)
if(clf.n == 2)
    plot_trial(2,clf,traj_u,traj_y) ;
elseif(clf.n == 3)
    plot_trial(3,clf,traj_u,traj_y) ;
else
    disp(['Skipping Trial #' num2str(j) ' #states: ' num2str(td.n_states) ' #o/p: ' num2str(td.n_outputs) ' #int: ' num2str(td.n_integrators)]) ;
end
end


function plot_trial(n_dim,clf,u,y)
figure ;
for i = 1:clf.n
    subplot(clf.n,1,i);
    plot(u(:,i)) ; grid on ; title(['u(' num2str(i) ')']);
end
    
    
figure ;
if(n_dim==2)
    plot(y(:,1), y(:,2), 'LineWidth', 5) ; hold on ;
    plot(clf.X0(1), clf.X0(2), 'k*') ;
else
    plot3(y(:,1), y(:,2), y(:,3), 'LineWidth', 5) ; hold on ;
    plot3(clf.X0(1), clf.X0(2), clf.X0(3), 'k*') ;
    zlabel('z') ;
end
xlabel('x'); ylabel('y') ; grid on ;
    
plot_boxes(n_dim, [clf.box_size(1:2:end,:); clf.box_size(2:2:end,:)],'red',0.3,false) ;
plot_boxes(n_dim, [clf.goal_size(1:2:end,:); clf.goal_size(2:2:end,:)],'green',0.3,true) ;
plot_boxes(n_dim, [clf.Y_min; clf.Y_max],'blue',0.01,false) ;
end
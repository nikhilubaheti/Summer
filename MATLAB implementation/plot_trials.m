function plot_trials(data)
    if(nargin < 1)
        data = loadjson('mydata.json') ;
    end
    close all ;
    
    data.trialconf
    for j=1:data.trialconf.number_trials
        td = data.trials{j} ;
        td.n_states = length(td.problem_instance.Xinit) ;
        td.n_outputs = size(td.problem_instance.U.H, 2) ;
        td.n_integrators = td.n_states / td.n_outputs ;
        td.trial_idx = j ;
        if(td.n_outputs == 2)
            plot_trial(2,td) ;
        elseif(td.n_outputs == 3)
            plot_trial(3,td) ;
        else
            disp(['Skipping Trial #' num2str(j) ' #states: ' num2str(td.n_states) ' #o/p: ' num2str(td.n_outputs) ' #int: ' num2str(td.n_integrators)]) ;
        end
    end
end


function plot_trial(n_dim,td)
    if(td.realizable)
        u_idx = 2 + [1:td.n_outputs] ;
        y_idx = u_idx(end) + [1:td.n_outputs] ;
        u = td.trajectory(:,u_idx) ;
        y = td.trajectory(:,y_idx) ;
        figure ; plot(u(:,1)) ; grid on ; title('u(1)') ;
    end
    
    
    figure ;
    if(n_dim==2)
        if(td.realizable)
            plot(y(:,1), y(:,2), 'LineWidth', 5) ; hold on ;
        end
        plot(td.problem_instance.Xinit(1), td.problem_instance.Xinit(2), 'k*') ;
    else
        if(td.realizable)
            plot3(y(:,1), y(:,2), y(:,3), 'LineWidth', 5) ; hold on ;
        end
        plot3(td.problem_instance.Xinit(1), td.problem_instance.Xinit(2), td.problem_instance.Xinit(3), 'k*') ;
        zlabel('z') ; 
    end
    xlabel('x'); ylabel('y') ; grid on ;
    title(['Trial #' num2str(td.trial_idx)]) ;
    
    plot_boxes(n_dim, td.problem_instance.obstacles,'red',0.3,false) ;
    plot_boxes(n_dim, td.problem_instance.goals,'green',0.3,true) ;
    plot_boxes(n_dim, {td.problem_instance.Y},'blue',0.01,false) ;
    
    % Check input and state
    if(td.realizable)
        n_obs = length(td.problem_instance.obstacles) ;
        input_violation = false ;
        output_violation = zeros(n_obs,1) ;
        for j=1:size(u,1)
            if(sum(td.problem_instance.U.H * u(j,:)' <= td.problem_instance.U.K') ~= 2*td.n_outputs)
                input_violation = true ;
                 break ;
            end
        end
        for j=1:size(u,1)
            for k=1:n_obs
                if(sum(td.problem_instance.obstacles{k}.H * y(j,:)' > td.problem_instance.obstacles{k}.K') == 2*td.n_outputs)
                    output_violation(k) = true ;
                    break ;
                end
            end
        end

        input_txt = 'Input passed' ;
        output_txt = 'Output passed' ;
        if(input_violation)
            input_txt = 'Input FAILED' ;
        end
        if(any(output_violation))
            output_txt = ['Output (' num2str(find(output_violation' == 1)) ') FAILED'] ;
        end

        disp(['Trial #' num2str(td.trial_idx) ' (n_int=' num2str(td.n_integrators) '): ' input_txt ',  ' output_txt]) ;
    else
        disp(['Trial #' num2str(td.trial_idx) ' (n_int=' num2str(td.n_integrators) '): Not realizable']) ;
    end
    
%     for j=1:length(td.problem_instance.obstacles)
%         V = con2vert(td.problem_instance.obstacles{j}.H, td.problem_instance.obstacles{j}.K') ;
%         min_V = min(V) ;
%         max_V = max(V) ;
%         if(n_dim==2)
%             h = prism_2D(min_V(1), min_V(2), max_V(1)-min_V(1), max_V(2)-min_V(2), 'red') ;
%         else
%             h = prism_3D(min_V(1), min_V(2), min_V(1), max_V(1)-min_V(1), max_V(2)-min_V(2), max_V(3)-min_V(3), 'red') ;
%         end
%     end
%     for j=1:length(td.problem_instance.goals)
%         V = con2vert(td.problem_instance.goals{j}.H, td.problem_instance.goals{j}.K') ;
%         min_V = min(V) ;
%         max_V = max(V) ;
%         h = prism_2D(min_V(1), min_V(2), max_V(1)-min_V(1), max_V(2)-min_V(2), 'green') ;
%         m = mean(V) ;
%         text(m(1),m(2), num2str(j)) ;
%     end
end

function plot_boxes(n_dim, objects,clr,alpha,txt_flag)
    for j=1:length(objects)
        V = con2vert(objects{j}.H, objects{j}.K') ;
        min_V = min(V) ;
        max_V = max(V) ;
        if(n_dim==2)
            h = prism_2D(min_V(1), min_V(2), max_V(1)-min_V(1), max_V(2)-min_V(2), clr,alpha) ;
        else
            h = prism_3D(min_V(1), min_V(2), min_V(1), max_V(1)-min_V(1), max_V(2)-min_V(2), max_V(3)-min_V(3), clr,alpha) ;
        end
        if(txt_flag)
            m = mean(V) ;
            if(n_dim==2)
                text(m(1),m(2), num2str(j)) ;
            else
                text(m(1),m(2),m(3),num2str(j),'HorizontalAlignment','left','FontSize',8);
            end
        end
    end
end


function h = prism_2D(x, y, w, l, clr,alpha)
    [X Y] = prism_faces_2D(x, y, w, l);

    faces(1, :) = [1:4];

    h = patch('Faces',faces,'Vertices',[X' Y'],'FaceColor',clr,'FaceAlpha',alpha) ;
end
function [X Y] = prism_faces_2D(x, y, w, l)
    X = [x x   x+w x+w];
    Y = [y y+l y+l y];
end


% Draw a 3D prism at (x, y, z) with width w,
% length l, and height h. Return a handle to
% the prism object.
function h = prism_3D(x, y, z, w, l, h, clr,alpha)
    [X Y Z] = prism_faces_3D(x, y, z, w, l, h);

    faces(1, :) = [4 2 1 3];
    faces(2, :) = [4 2 1 3] + 4;
    faces(3, :) = [4 2 6 8];
    faces(4, :) = [4 2 6 8] - 1;
    faces(5, :) = [1 2 6 5];
    faces(6, :) = [1 2 6 5] + 2;

    h = patch('Faces',faces,'Vertices',[X' Y' Z'],'FaceColor',clr,'FaceAlpha',alpha) ;
end



% Compute the points on the edge of a prism at
% location (x, y, z) with width w, length l, and height h.
function [X Y Z] = prism_faces_3D(x, y, z, w, l, h)
    X = [x x x x x+w x+w x+w x+w];
    Y = [y y y+l y+l y y y+l y+l];
    Z = [z z+h z z+h z z+h z z+h];
end

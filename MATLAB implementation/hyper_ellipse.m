classdef hyper_ellipse
    properties
        center = [];
        axis = [];
        const = 0;
        degree = 2;
    end
    methods
        function ellipse = ellips_obs_constr(h_ellipse,clf,box_size, goal, n, obs_flag,X0)
            % '''
            % Generate an ellipsoid-shaped obstacle to approximate a box
            % box_size: a list containing the minimum and maximum of each dimension
            % goal:     a list containing each goal position
            % n:    the dimension of the position space
            % obs_flag:	the flag determine whether the box should be inscribed(1) or circumscribed(0)
            % '''
            ellipse = h_ellipse;
            goal_num = size(goal,2);
            goal(:,goal_num+1) = X0(1:n);
            goal_num = goal_num+1;
            if obs_flag == 0
                p = 1.0;
            else
                p = 1.0;               ...  % the initial degree of this super ellipsoid
            end
        if clf.commands.fixed_p 
            p = 16;
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
                if n == 2
                    plot_se(box_center,dist*(n*clf.commands.ellipse_const_mul^(1/max_p)),max_p); hold on;
                    plot_boxes(n, [box_size(1:2:end,:); box_size(2:2:end,:)],'red',0.3,true) ;
                    plot(goal(1,:),goal(2,:),'k*');hold off;
                end
                return;
                end
                a = dist.^p;
                flag = true;
                
                for j = 1:goal_num
                    exo = box_center - goal(:,j);
                    g = (exo.^p)./a;
                    g = sum(g);
                    if g <= n*clf.commands.ellipse_const_mul && obs_flag == 0
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
%         scale = min(a);
%         a = a./scale;
        if obs_flag == 0
            b = n*clf.commands.ellipse_const_mul;
        else
            b = 1;%scale;
        end
        
        ellipse.center = box_center;
        ellipse.axis = a;
        ellipse.const = b;
        ellipse.degree = p;
        end
    end
end
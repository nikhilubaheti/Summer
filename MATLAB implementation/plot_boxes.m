function plot_boxes(n_dim, objects,clr,alpha,txt_flag)
    for j=1:size(objects,2)
        min_V = objects(1:n_dim,j) ;
        max_V = objects(n_dim+1:end,j) ;
        if(n_dim==2)
            h = prism_2D(min_V(1), min_V(2), max_V(1)-min_V(1), max_V(2)-min_V(2), clr,alpha) ;
        else
            h = prism_3D(min_V(1), min_V(2), min_V(3), max_V(1)-min_V(1), max_V(2)-min_V(2), max_V(3)-min_V(3), clr,alpha) ;
        end
        if(txt_flag)
            m = (min_V+max_V)/2 ;
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

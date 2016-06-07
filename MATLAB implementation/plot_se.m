% Plots Superellipse  (actually only works for 2D for now)
% Inputs:
%   C = [xc;yc]
%   R = [rx;ry]
%   n = order
% Details:
%   Supperelippse equation:  ( (x-xc)/rx )^n + ( (y-yc)/ry )^n = 1
%   If rx=ry=r, and n=2, then get circle with radius r, centered at (xc,yc)
%   If n=2, get ellipse with x and y axis as rx, ry.
%   If n=4, get rectellipse
function plot_se(C,R,n)
    ang=0:0.01:2*pi; 
    
    xp=R(1)*cos(ang).^(2/n);
    yp=R(2)*sin(ang).^(2/n);
    
    plot(C(1)+xp,C(2)+yp);
    grid on; axis equal ;
end
% Plots Superellipse  (actually only works for n=2 for now)
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
    ang=0:0.5:360; 
    xp=R(1)*abs(cosd(ang)).^(2/n).*sign(cosd(ang));
    yp=R(2)*abs(sind(ang)).^(2/n).*sign(sind(ang));
%     figure;
    plot(C(1)+xp,C(2)+yp);
%     grid on; axis equal ;
end
% Superellipse  (Use plot_se() to visualize)
% Inputs:
%   p = [x;y]
%   dp = [dx;dy]
%   Lf2p, LgLfp s.t. d2p = Lf2p + LgLfp*u
%   C = [xc;yc]
%   R = [rx;ry]
%   n = order
%   sgn = -1 for reversing inequality
% Details:
%   Supperelippse equation:  ( (x-xc)/rx )^n + ( (y-yc)/ry )^n = 1
%   If rx=ry=r, and n=2, then get circle with radius r, centered at (xc,yc)
%   If n=2, get ellipse with x and y axis as rx, ry.
%   If n=4, get rectellipse
% Outputs:
%   g = ( (x-xc)/rx )^n + ( (y-yc)/ry )^n - 1  (g >= 0)
%   dg = time  derivative of g
%   Lf2_g LgLfg s.t. d2g = Lf2_g + LgLf_g*u
function [g dg Lf2_g LgLf_g] = calc_se(p,dp,Lf2p,LgLfp,C,R,n, sgn)
    pCR = (p-C)./R ;
    g = sum( pCR.^n ) - 1 ;
    dg = ( n*pCR.^(n-1) )' * dp ;
    Lf2_g = ( n*(n-1)*(pCR).^(n-2) )' * (dp.*dp) + ...
            ( n*pCR.^(n-1) )' * Lf2p ;
    LgLf_g = ( n*pCR.^(n-1) )' * LgLfp ;
    
    % flip signs if needed
    if(sgn == -1)
        g = -g ;
        dg = -dg ;
        Lf2_g = -Lf2_g ;
        LgLf_g = -LgLf_g ;
    end
end

function Robust_CBF_CLF_Pointmass
clear;clc;close all;
m=1;
degree=2;   % degree where u is applied
amp=1;
w=1/10;
amp = 0.5;
w = 1;
k=10;
x0=0.3; y0=0.5; z0=2; r0=0.5; delta=ones(3,1)/m;
F = [zeros(3*(degree-1),3),eye(3*(degree-1));zeros(3,3*degree)];
G = [zeros(3*(degree-1),3);eye(3)/m];
[K,P,E] = lqr(F,G,eye(size(F,2)),eye(size(G,2)));
x_init = traj_gen(0);
x_init = x_init + 0.5*x_init;%.*rand(size(x_init));
%opts1 = odeset('RelTol',1e-5,'AbsTol',1e-6);
Q = -((F-G*K)'*P+P*(F-G*K));
[t1, x1] = ode45(@continous_dynamics, [0:0.1:50], x_init);
save('x1.mat','x1','t1');
plots();

    function plots()
        figure(1);  
        x_d = arrayfun(@traj_gen, t1, 'UniformOutput', false);
        x_d = cell2mat(x_d);
        x_d = reshape(x_d,3*degree,numel(x_d)/(3*degree))';

        plot(t1,x1(:,1:3),'*')
        hold on;
        plot(t1,x_d(:,1:3));
        hold off;
        title('x1-xd');
        figure(2);
        plot(t1,x_d(:,1:3)-x1(:,1:3));
        title('Error');
        legend('x','y','z');
        figure;
        animate();
        title('Robust vs Non-Robust CLF-CBF Controller');
    end

    function [x_dot] = continous_dynamics(t,x)
        xd = traj_gen(t);
        %u = CLF_controller(t,xd,x);
        u = exp_CLFCBF_controller(t,xd,x);
        %u = linear_controller(t,xd,x);
        noise = 0.9*sin(t)*ones(3,1);        
        disp(t);
        x_dot = F*x+G*(u + 0*noise);
    end

    function u = linear_controller(t,xd,x)
        % FF term only for circular trajectory
        [u_ff, xd_dot] = feedforward(t,xd);
        u = K*(xd-x) + u_ff;
    end

    function u = CLFCBF_controller(t,xd,x)
%         [u_ff, xd_dot] = feedforward(t,xd);
%         u_ff =  - m*[xd(1); xd(2); 0];
%         V = (x-xd)'*P*(x-xd);
%         dV_dx = 2*(x-xd)'*P;
%         LfV = dV_dx*F*(x);
%         LgV = dV_dx*G;
%         gamma=1;
%         b = -LfV - LgV*u_ff + dV_dx*xd_dot -(x-xd)'*Q*(x-xd); % - gamma*V; %
%         A = [LgV -1];
%         [B,dB_dx] = fcn_CBF_PointMass(x,[x0 y0 z0 r0]);
%         LfB = dB_dx*F*x;
%         LgB = dB_dx*G;
%         A = [A; LgB 0];
%         b = [b; -LfB - LgB*u_ff + (gamma*B)];
%         W=diag([1,1,1,1000]);
%         %opts = optimset('Algorithm','interior-point-convex','Display','off');
%         opts = optimset('LargeScale','off','Display','off');
%         Aeq = []; beq=[];
%         opt = quadprog(W,[0;0;0;0],A,b,Aeq,beq,[],[],[],opts);
%         u = opt(1:3) + u_ff;
% %         H = (x(1)-x0)^2 + (x(2)-y0)^2 + (x(3)-z0)^2 - r0^2;
% %         figure(5);
% %         plot(t,V,'r*');
% %         title('V');
% %         hold on;
% %         figure(6);
% %         plot(t,B,'b*')
% %         hold on;
% %         title('B');
% %         figure(7);
% %         plot(t,norm(u),'g*')
% %         hold on;
% %         title('Control Input');
    end


function u = CLFCBF_controller1(t,xd,x)
        alpha=1; e=1;
        FF =  - m*[xd(1); xd(2); 0];
        V_weight = [alpha/2     0     0     e/2  0   0;
                       0    alpha/2   0      0  e/2  0;
                       0       0    alpha/2  0   0  e/2;
                       e/2     0      0     m/2  0   0;
                       0       e/2    0      0  m/2  0;
                       0       0      e/2    0   0  m/2;];       
        V = (x-xd)'*V_weight*(x-xd);
        dV_dx = 2*V_weight*(x-xd);
        x_dot = [x(4:6,1); 1/m]; % Missing u
        gamma = 1;
        b = -dV_dx(1:3)'*x_dot(1:3) - gamma*V - dV_dx(4:6)'*FF/m;
        A = [dV_dx(4:6)'/m, -1];
        
        g = (x(1)-x0)^2 + (x(2)-y0)^2 + (x(3)-z0)^2 - r0^2;
        dg_dx = [2*(x(1)-x0); 2*(x(2)-y0); 2*(x(3)-z0); 0; 0; 0];
        h = g + dg_dx(1:3)'*x(4:6);
        B = 1/h;
        dB_dx =(-1/(h^2))* [2*(x(1)-x0)+2*x(4); 2*(x(2)-y0)+2*x(5); 2*(x(3)-z0)+2*x(6); ...
                                   2*(x(1)-x0); 2*(x(2)-y0); 2*(x(3)-z0)];
        gamma_B = 10;
        A = [A; [dB_dx(4:6)' 0]];
% TODO check 1/m;        
        b = [b; -dB_dx(1:3)'*x(4:6) + (gamma_B/B) - dB_dx(4:6)'*FF/m];
        W=diag([1,1,1,1000]);
        
        %opts = optimset('Algorithm','interior-point-convex','Display','off');
        opts = optimset('LargeScale','off','Display','off');
        opt = quadprog(W,[0;0;0;0;],A,b,[],[],[],[],[],opts);
        u = opt(1:3) + FF;
        
%        figure(1);
%        plot(t,u(1),'r*');
%        hold on;
end

    function u = exp_CLFCBF_controller(t,xd,x)
%        [u_ff, xd_dot] = feedforward(t,xd);
%         u_ff =  - m*[xd(1); xd(2); 0];
%         V = (x-xd)'*P*(x-xd);
%         dV_dx = 2*(x-xd)'*P;
%         LfV = dV_dx*F*(x);
%         LgV = dV_dx*G;
%         gamma=1;
%         b = -LfV + dV_dx*xd_dot  - LgV*u_ff - (x-xd)'*Q*(x-xd); % - gamma*V;
%         A = [LgV -1];
%         [B,dB_dx] = fcn_expCBF_PointMass(x,[x0 y0 z0 r0]);
%         LfB = dB_dx*F*x;
%         LgB = dB_dx*G;
%         A = [A; -LgB 0];
%         b = [b; LfB + LgB*u_ff + (gamma*B)];
%         W=diag([1,1,1,100]);
%         %opts = optimset('Algorithm','interior-point-convex','Display','off');
%         opts = optimset('LargeScale','off','Display','off','TolFun',1e-9);
%         Aeq = []; beq=[];
%         opt = quadprog(W,[0;0;0;0],A,b,Aeq,beq,[],[],[],opts);
%         u = opt(1:3) + u_ff;
        [u_ff, xd_dot] = feedforward(t,xd);
        V = (x-xd)'*P*(x-xd);
        dV_dx = 2*(x-xd)'*P;
        LfV = dV_dx*F*(x);
        LgV = dV_dx*G;
        gamma=10;
        b = -LfV + dV_dx*xd_dot  - LgV*u_ff - (x-xd)'*Q*(x-xd); % - gamma*V;
        A = [LgV -1];
        [B,dB_dx] = fcn_expCBF_PointMass(x,[x0 y0 z0 r0]);
        LfB = dB_dx*F*x;
        LgB = dB_dx*G;
        A = [A; -LgB 0];
        b = [b; LfB + LgB*u_ff + (gamma*B)];
        W=diag([1,1,1,1000]);
        %opts = optimset('Algorithm','interior-point-convex','Display','off');
        opts = optimset('LargeScale','off','Display','off','TolFun',1e-9);
        Aeq = []; beq=[];
        opt = quadprog(W,[0;0;0;0],A,b,Aeq,beq,[],[],[],opts);
        u = opt(1:3) + u_ff;
    end

    function u = CLF_controller(t,xd,x)
        [u_ff, xd_dot] = feedforward(t,xd);
        V = (x-xd)'*P*(x-xd);
        dV_dx = 2*(x-xd)'*P;
        LfV = dV_dx*F*(x);
        LgV = dV_dx*G;
        gamma=12;
        b = -LfV + dV_dx*xd_dot  - LgV*u_ff - (x-xd)'*Q*(x-xd); - gamma*V;
        %b = -LfV - gamma*V + LgV*u_ff;
        A = LgV;
        Aeq=[];
        beq=[];
        opts = optimset('LargeScale','off','Display','off','TolFun',1e-6);
%         figure(5);
%         plot(t,u(1:3),'*');
%         hold on;
        figure(6);
        plot(t,V,'r*');
        hold on;
                    %u = quadprog(eye(3),[0;0;0],Aeq,beq,A,b,[],[],[0;0;0],opts) + u_ff;
%        if norm(LgV)>0.001
%            u = quadprog(eye(3),[0;0;0],Aeq,beq,A,b,[],[],[0;0;0],opts) + u_ff;
%        else 
%            u = u_ff;
%        end
    end

    function x_d = traj_gen(t)
        % assuming 3-d
        x_d_full = [ amp*sin(t*w), amp*cos(t*w), t/k, amp*w*cos(t*w), -amp*w*sin(t*w), 1/k, -amp*w^2*sin(t*w), -amp*w^2*cos(t*w), 0, -amp*w^3*cos(t*w), amp*w^3*sin(t*w), 0, amp*w^4*sin(t*w), amp*w^4*cos(t*w), 0, amp*w^5*cos(t*w), -amp*w^5*sin(t*w), 0, -amp*w^6*sin(t*w), -amp*w^6*cos(t*w), 0, -amp*w^7*cos(t*w), amp*w^7*sin(t*w), 0, amp*w^8*sin(t*w), amp*w^8*cos(t*w), 0, amp*w^9*cos(t*w), -amp*w^9*sin(t*w), 0];
        x_d = x_d_full(1:3*degree)';
    end
    
    function [u_ff, xd_dot] = feedforward(t,xd)
      x_d_full = [ amp*sin(t*w), amp*cos(t*w), t/k, amp*w*cos(t*w), -amp*w*sin(t*w), 1/k, -amp*w^2*sin(t*w), -amp*w^2*cos(t*w), 0, -amp*w^3*cos(t*w), amp*w^3*sin(t*w), 0, amp*w^4*sin(t*w), amp*w^4*cos(t*w), 0, amp*w^5*cos(t*w), -amp*w^5*sin(t*w), 0, -amp*w^6*sin(t*w), -amp*w^6*cos(t*w), 0, -amp*w^7*cos(t*w), amp*w^7*sin(t*w), 0, amp*w^8*sin(t*w), amp*w^8*cos(t*w), 0, amp*w^9*cos(t*w), -amp*w^9*sin(t*w), 0];
      u_ff = m*x_d_full(3*degree+1:3*degree+3)';
      xd_dot = F*xd + G*u_ff;
    end

    function animate
        plot3(x1(:,1),x1(:,2),x1(:,3));
        ax=1;
        axis([-ax ax -ax ax 0 10]);
        hold on;
        [xs,ys,zs] = sphere;
        h = surfl(xs*r0+x0,ys*r0+y0,zs*r0+z0);
        set(h, 'FaceAlpha', 0.5)
        shading interp      
        hold off;
    end
end
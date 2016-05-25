classdef CLF_CBF_QP
    %      '''
    %         Initialize the CLF_CBF problem formulation
    %         intopic: where the state is read from
    %         Outtopic: where the input must be published to
    %         m: integrator chain order
    %         n: dimension of output space = number of control inputs
    %         A(m*n,m*n): State matrix for the given integrator chain
    %         B(m*n,n): Control matrix for the given integrator chain
    %         target: the desired goal location where the state must be driven to.
    %         	This is updated to the next goal location once this state converges to the target
    %         box_size: Contains the extremums of all the obstacles in each dimension
    % 		U_min: Contains the minimum value of the control input for each dimension
    % 		U_max: Contains the maximum value of the control input for each dimension
    % 		Y_min: Contains the minimum value of the position for each dimension
    % 		Y_max: Contains the maximum value of the position input for each dimension
    % 		X0(Currently unused): The initial state vector
    % 		cbf_prms: Contains B, dB,..., d^(m-1)B, Lf^(m)B,LgLf^(m-1)B terms in symbolic form for each obstacle
    % 		Ks: contains the feedback matrices for each obstacle
    % 		Realizability: Set to False if there problem may not be solved for some reason, terminates the problem
    % 		cbf_Y: Contains B, dB,..., d^(m-1)B, Lf^(m)B,LgLf^(m-1)B terms in symbolic form for the outter bounds on the output vector
    % 		K_Ys: Contains the feedback matrices for the boundary box
    % 		'''
    properties (Access = public)
        %         Initializations
%         error = None
        m = 2;
        n = 2;
        A = [];
        B = [];
        U_min = [];
        U_max = [];
        Y_min = [];
        Y_max = [];
        Realizability = true;
        box_size = [];
        obs_num = 0;
        Ks = [];
        cbf_p = [];
        K_Ys = [];
        goal_size = [];
        Reached_Goal = false;
        goal_num = 0;
        X0 = [];
        cbf_Y = [];
        Q = [];
        R = [];
        target_state = [];
        T = [];
        vel = [];
        P = [];
        goal = [];
        u = [];
    end
    
    methods
        function clf = initialize(clf)
            clf.A = diag(ones(1,(clf.m-1)*clf.n),clf.n);
            clf.B = [zeros((clf.m-1)*clf.n, clf.n);eye(clf.n)];
            clf.U_min = -ones(clf.n);
            clf.U_max = ones(clf.n);
            clf.Y_min = -10*ones(clf.n);
            clf.Y_max = 10*ones(clf.n);
            clf.X0 = zeros(clf.n*clf.m);
            clf.Q = eye(clf.m*clf.n);
            clf.R = 10*eye(clf.n);
            clf.target_state = zeros(clf.m*clf.n);
            clf.T = min(4,25/clf.goal_num);
            clf.vel = zeros(clf.n,1);
            clf.goal = zeros(clf.n,1);
            [~,clf.P,~] = lqr(clf.A,clf.B,clf.Q,clf.R);
        end
        
        function u = controller(clf,t,x)
            gamma = 1;
            [A_ineq_lyap,b_ineq_lyap] = clf.lyap_constraints(t,x,gamma);
            [A_ineq_ip,b_ineq_ip] = clf.control_ip_constraints();
            [A_ineq_bound,b_ineq_bound] = clf.state_bound_constraints(x);
            [A_ineq_obs,b_ineq_obs] = clf.obs_constraints(x);
            A_ineq = [A_ineq_lyap;A_ineq_ip;A_ineq_bound;A_ineq_obs];
            b_ineq = [b_ineq_lyap;b_ineq_ip;b_ineq_bound;b_ineq_obs];
%             Normalize every inequality constraint
            norm = sqrt(sum(abs(A_ineq).^2,2));
            A_ineq = A_ineq./repmat(norm,1,clf.n+1);
            b_ineq = b_ineq./norm;
%           '''
%           Quadratic cost matrix: diag([ones(n),100])
%           high cost to slack variable to ensure we do not collide with obstacles
%           '''
            C_qp = eye(clf.n+1);
            C_qp(clf.n+1,clf.n+1) = 100;
            options = optimoptions('quadprog','Display','off');
            u_slack = quadprog(C_qp,zeros(clf.n+1,1),A_ineq,b_ineq,[],[],[],[],[],options);
            if length(u_slack) ~= clf.n+1
                fprintf('QP failed\n');
                u = zeros(clf.n,1);
            else
                u = u_slack(1:clf.n);
                u = max(clf.U_min,min(u,clf.U_max));
            end
    end
    
    function [A_ineq,b_ineq] = lyap_constraints(clf,t,x,gamma)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %         Input:
    %             clf: clf_cbf object
    %             t: current simulation time
    %             x(m*n,1): state vector
    %         Output:
    %             A_ineq(1,n+1): Inequality matrix A for the input u
    %             b_ineq(1): Inequality constant
    %             A_ineq*[u;delta] <= b_ineq
    %         Tunning variable:
    %             gamma: convergence rate for exponential stability
    %
    %         V = x'*P*x
    %         F: state vector = A*X
    %         G: Control vector = B
    %         '''
    %         CLF_QP:LgV*u + delta<= -gamma*V-LfV
    %         LfV = (dV/dx)*F*x
    %         LgV = (dV/dx)*G
    %         delta: slack variable
    %         '''
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    F = clf.A;
    G = clf.B;
%     clf.state_vec = x;
    pos_d = clf.X0(1:clf.n) + clf.vel*min(t,clf.T);
    xd = [pos_d; zeros((clf.m-1)*clf.n,1)];
    clf.Reached_Goal = chk_collision(x(1:clf.n),clf.goal_size);
    if clf.Reached_Goal
        fprintf('reached goal!!');
    end
    V = (x-xd)'*clf.P*(x-xd);
    dV_dx = 2*(x-xd)'*clf.P;
    LfV = dV_dx*F*x;
    LgV = dV_dx*G;
    b_ineq = -LfV - gamma*V;
    A_ineq = [LgV,1];
    end
    
    function [A_ineq,b_ineq] = control_ip_constraints(clf)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %         Input:
    %             clf: clf_cbf object must contain U_min and U_max values
    %             t: current simulation time
    %             x(m*n,1): state vector
    %         Output:
    %             A_ineq(2*n,n+1): Inequality matrix A for the input u
    %             b_ineq(2*n): Inequality constant
    %             A_ineq*[u;delta] <= b_ineq
    %
    %         '''
    %         U_QP: Contraints on control input U
    %         I*u >= U_min
    %         I*u <= U_max
    %         '''
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    A_ineq = [-eye(clf.n), zeros(clf.n,1);eye(clf.n), zeros(clf.n,1)];
    b_ineq = [-clf.U_min;clf.U_max];
    if size(b_ineq,1) ~= 2*clf.n;
        fprintf('b_ineq not of correct size for input(u) inequality constraints');
    end
    end
    
    function [A_ineq,b_ineq] = state_bound_constraints(clf,x)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %         Input:
    %             clf: clf_cbf object must contain feedback gains K_min and
    %             K_max
    %             t: current simulation time
    %             x(m*n,1): state vector
    %         Output:
    %             A_ineq(2*n,n+1): Inequality matrix A for the input u
    %             b_ineq(2*n): Inequality constant
    %             A_ineq*[u;delta] <= b_ineq
    %        Virtual Tunning parameter:
    %             K: feedback matrix for eta vector; actually found using pole
    %             placement
    %
    %         '''
    %         CBF_QP:  - LgLf^(m-1)B*u <= Lf^(m)B - K*eta
    %         eta = [B(x);dB(x);...;d^(m-1)B(x)];
    %         B(x) = sum [(xi/ai)^p] - 1
    %         Constraint to be outside obstacles modelled as CBF
    %         '''
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    A_ineq_min = [-eye(clf.n),zeros(clf.n,1)];
    A_ineq_max = [eye(clf.n),zeros(clf.n,1)];
    eta_min = x - [clf.Y_min; zeros(clf.n*(clf.m-1),1)];
    b_ineq_min = sum(clf.K_Ys(1:2:end,:).*reshape(eta_min,clf.n,clf.m),2);
    eta_max = x - [clf.Y_max; zeros(clf.n*(clf.m-1),1)];
    b_ineq_max = -sum(clf.K_Ys(2:2:end,:).*reshape(eta_max,clf.n,clf.m),2);
    A_ineq = [A_ineq_min;A_ineq_max];
    b_ineq = [b_ineq_min;b_ineq_max];
    end
    
    function [A_ineq,b_ineq] = obs_constraints(clf,x)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %         Input:
    %             clf: clf_cbf object must contain feedback gains and obstacle
    %             eta vectors
    %             t: current simulation time
    %             x(m*n,1): state vector
    %         Output:
    %             A_ineq(2*n,n+1): Inequality matrix A for the input u
    %             b_ineq(2*n): Inequality constant
    %             A_ineq*[u;delta] <= b_ineq
    %        Virtual Tunning parameter:
    %             K: feedback matrix for eta vector; actually found using pole
    %             placement
    %
    %         '''
    %         Y_QP: Contraints on state position vector outter bounds
    %         For minimum bounds: x >= x_min
    %         B(x) = x-x_min
    %         eta = [B(x);dB(x);...;d^(m-1)B(x)];
    %         For integrator chain case:
    %           -u <= K*eta
    %
    %
    %         For maximum bounds: x <= x_max
    %         B(x) = x-x_max
    %         eta = [B(x);dB(x);...;d^(m-1)B(x)];
    %         For integrator chain case:
    %           u <= -K*eta
    %         '''
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    b_ineq = [];
    A_ineq = [];
    for i = 1:clf.obs_num
%         cbf = clf.cbf_p(i);
        cbf_sub = clf.cbf_p(i).subs_cbf(x);
        if cbf_sub.eta(1)<0
            fprintf('in collision with obstacle');
        end
        K = clf.Ks(i,:);
        b_ineq_cbf = K*cbf_sub.eta+cbf_sub.LfB;
        if length(b_ineq_cbf) ~= 1
            fprintf('b_ineq is not a value');
        end
        b_ineq = [b_ineq;b_ineq_cbf];
        A_ineq = [A_ineq;-cbf_sub.LgB,0];
    end
    end

    end
end
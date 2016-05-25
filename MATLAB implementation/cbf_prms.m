classdef cbf_prms
    properties
        eta = [];
        LfB = [];
        LgB = [];
        state = [];
    end
    methods
        function cbf = cbf_gen(cbf_p,m, n, geom)
            % '''
            % Generate a CBF and its higher order derivative up to nth order
            % including Lf^n B and LgLf^(n-1)B
            % Return the list containing B, dB, ..., d^n B, Lf^n B, Lg Lf^(n-1 B and the
            %     state name
            % All the types are converted to Matrix in sympy
            %
            % m       the order of the chain of integrator
            % n       the dimension of position space
            % obs_type    a string indicating the obstacle's type
            % geom        a dictionary storing the geometric parameters of this obstacle
            %
            % for sphere obstacle
            %     ------------------------------------------------------------------------------------
            %     geom['center']  a list containing the position of the center
            %     geom['radius']  the radius of this sphere
            %
            %     0
            %     ------------------------------------------------------------------------------------
            %     '''
            fprintf('Generate expression for CBF using Sympy...\n');
            statex = [sym('x%d',[n,1]);reshape(sym('d%dx%d',[m-1,n],'real')',(m-1)*n,1)];
            func = sum(((sym('x%d',[n,1])-geom.center)./geom.axis).^geom.degree)-geom.const;
            cbf = cbf_p;   % create a list for CBF
            cbf.eta = sym('eta',[m,1]);
            % cbf.eta = zeros(m,1);
            cbf.eta(1) = func;
            % dstate = reshape(sym('d%dx%d',[m,n],'real')',m*n,1);
            dstate = [statex(n+1:end);zeros(n,1)];
            for i = 2:m
                cbf.eta(i) = jacobian(cbf.eta(i-1),statex)*dstate;
            end
            % compute the lie derivative of the highest derivatives
            cbf.LfB = jacobian(cbf.eta(m),statex)*dstate;
            cbf.LgB = jacobian(cbf.eta(m),statex((m-1)*n+1:end));
            cbf.state = statex;
        end
        
        
        function cbf_sub = subs_cbf(cbf_p, x)
            % """
            % Substitute the value of current state and return
            % the numerical value as a list
            %
            % cbf_prms    the expressions for B and its higher derivatives
            % statec      a list or tuple containing the values of current state
            % """
            cbf_sub = cbf_p;
            % create a dictionary for substitution
            cbf_sub.eta = double(subs(cbf_p.eta,cbf_p.state,x));
            cbf_sub.LfB = double(subs(cbf_p.LfB,cbf_p.state,x));
            cbf_sub.LgB = double(subs(cbf_p.LgB,cbf_p.state,x));
            %             '''
            %             structure of cbf_val: is a list similar to cbf_prms, but only contain n+2 entries:
            %             each entry is just a floating point number
            %             entry 1-n: B, dB, d2B, ..., d(n-1)B, each one of them is a scalar
            %                 entry n+1: Lf^n B , the scalar term in dnB inpendent of u
            %             entry n+2: LgLf^(n-1)B, the 3 dim vector in dnB which is the coefficient
            %             of u in QP
            %             The structure of cbf_val is a bit messy at the moment, I will make it better later
            %             with better proficiency in python
            %             '''
        end
    end
end
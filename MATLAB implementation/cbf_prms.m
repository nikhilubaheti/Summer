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
            % including Lf^m B and LgLf^(m-1)B
            % Return the list containing B, dB, ..., d^m B, Lf^m B, Lg Lf^(m-1) B and the
            %     state name
            % 
            % m       the order of the chain of integrator
            % n       the dimension of position space
            % geom    ellipse geometry that contains the ellipse parameters
            %
            fprintf('Generate expression for CBF using Sympy...\n');
            statex = [sym('x%d',[n,1]);reshape(sym('d%dx%d',[m-1,n],'real')',(m-1)*n,1)];
            func = sum(((sym('x%d',[n,1])-geom.center)./geom.axis).^geom.degree)-geom.const;
            cbf = cbf_p;   % create a list for CBF
            cbf.eta = sym('eta',[m,1]);
            cbf.eta(1) = func;
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
            % cbf_p       the expressions for B and its higher derivatives
            % x           current state
            % """
            cbf_sub = cbf_p;
            cbf_sub.eta = double(subs(cbf_p.eta,cbf_p.state,x));
            cbf_sub.LfB = double(subs(cbf_p.LfB,cbf_p.state,x));
            cbf_sub.LgB = double(subs(cbf_p.LgB,cbf_p.state,x));
        end
    end
end
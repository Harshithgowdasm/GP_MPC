classdef VanDerPol_dynamics
    properties
        delta_t
        mu
        alpha
        beta
    end
    
    methods
        function obj = VanDerPol_dynamics(delta_t, mu, alpha, beta)
                beta = 0.8;
                alpha = 5;
                mu = 2;
                delta_t = 0.2;
            obj.delta_t = delta_t;
            obj.mu = mu;
            obj.alpha = alpha;
            obj.beta = beta;
        end
        
        function dx = f(obj, x)
            dx = [2*x(2); obj.mu*(1-obj.alpha*x(1)*x(1))*x(2)-obj.beta*x(1)];
        end
        
        function dx = f_u(obj, x, u)
            dx = obj.f(x) + [0; u];
        end
        
        function k1 = k1(obj, x, u)
            k1 = obj.f_u(x, u);
        end
        
        function k2 = k2(obj, x, u)
            k2 = obj.f_u(x + obj.k1(x, u) * obj.delta_t / 2, u);
        end
        
        function k3 = k3(obj, x, u)
            k3 = obj.f_u(x + obj.k2(x, u) * obj.delta_t / 2, u);
        end
        
        function k4 = k4(obj, x, u)
            k4 = obj.f_u(x + obj.k1(x, u) * obj.delta_t, u);
        end
        
        function x_next = f_ud(obj, x, u)
            x_next = x + (obj.delta_t / 6) * (obj.k1(x, u) + 2 * obj.k2(x, u) + 2 * obj.k3(x, u) + obj.k4(x, u));
        end
    end
end


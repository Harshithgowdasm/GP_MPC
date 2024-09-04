classdef VanDerPol_uncertanity_dynamics
    properties
        delta_t
        mu
        alpha
        beta
        gamma
    end
    % +0.001*randn(1)
    methods
        function obj = VanDerPol_uncertanity_dynamics(delta_t, mu, alpha, beta,gamma)
                beta = 0.8;%+0.8*0.05*randn(1);
                alpha = 5; %+5*0.05*randn(1);
                mu = 2; %+2*0.05*randn(1);
                gamma = 2; %+2*0.01*randn(1);
                delta_t = 0.2;
            obj.delta_t = delta_t;
            obj.mu = mu;
            obj.alpha = alpha;
            obj.beta = beta;
            obj.gamma = gamma;
        end
        
        function dx = f(obj, x,noise)
            dx = [(obj.gamma+2*0.01*noise)*x(2); (obj.mu+2*0.01*noise)*(1-(obj.alpha+5*0.01*noise)*x(1)*x(1))*x(2)-(obj.beta+0.8*0.05*noise)*x(1)];
        end
        
        function dx = f_u(obj, x, u,noise)
            dx = obj.f(x,noise) + [0; u];
        end
        
        function k1 = k1(obj, x, u,noise)
            k1 = obj.f_u(x, u,noise);
        end
        
        function k2 = k2(obj, x, u,noise)
            k2 = obj.f_u(x + obj.k1(x, u,noise) * obj.delta_t / 2, u,noise);
        end
        
        function k3 = k3(obj, x, u,noise)
            k3 = obj.f_u(x + obj.k2(x, u,noise) * obj.delta_t / 2, u,noise);
        end
        
        function k4 = k4(obj, x, u,noise)
            k4 = obj.f_u(x + obj.k1(x, u,noise) * obj.delta_t, u,noise);
        end
        
        function x_next = f_ud(obj, x, u,noise)
            x_next = x + (obj.delta_t / 6) * (obj.k1(x, u,noise) + 2 * obj.k2(x, u,noise) + 2 * obj.k3(x, u,noise) + obj.k4(x, u,noise));
        end
    end
end


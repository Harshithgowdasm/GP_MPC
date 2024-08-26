function [K] = KernalCov(x,z,hyp)

% n = 20;
% disp('x = gpml_randn(0.3, n, 1);')
% x =gpml_randn(0.3, n, 1);
% ell = 0.4;
%  sf = 1; 
%  z = 2;

% hyp= log([ell; sf]);
covfunc = {@covSEard};

% disp('K = feval(covfunc{:}, hyp.cov, x);')
K = feval(covfunc{:},hyp,x,z);

end


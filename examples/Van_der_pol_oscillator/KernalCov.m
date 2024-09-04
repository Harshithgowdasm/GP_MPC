function [K] = KernalCov(x,z,hyp)

% hyp= log([ell; sf]);
covfunc = {@covSEard};

% disp('K = feval(covfunc{:}, hyp.cov, x);')
K = feval(covfunc{:},hyp,x,z);

end


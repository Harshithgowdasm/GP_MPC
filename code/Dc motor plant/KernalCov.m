function [K] = KernalCov(x,z,hyp)


covfunc = {@covSEard};
K = feval(covfunc{:},hyp,x,z);

end


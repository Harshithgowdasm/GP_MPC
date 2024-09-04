%---------------------------------------------------------------------------------------------
% For Paper, 
% "A Safe Bayesian Optimization Algorithm for Tuning the Optical Synchronization System at European XFEL"
% by Jannis O. Lübsen, Maximilian Schütte, Sebastian Schulz, Annika Eichler
% Copyright (c) Institute of Control Systems, Hamburg University of Technology. All rights reserved.
% Licensed under the GPLv3. See LICENSE in the project root for license information.
% Author(s): Jannis Lübsen
%--------------------------------------------------------------------------------------------

function [hyp_new,nlZ,dnlZ]=fitGP(hyp, inf_, mean_, cov_, lik_, x, y, opts, varargin)
    % change options in the "opt" struct. E.g. opt.rndmInt.mean=...
    oldOpts.mode = 2;
    oldOpts.repeat = 10;
    oldOpts.rndmInt.mean = repmat([-20,50],size(hyp.mean,1),1);
    oldOpts.rndmInt.cov = repmat([1e-3,5],size(hyp.cov,1),1);
    oldOpts.rndmInt.lik = repmat([1e-3,1],size(hyp.lik,1),1);

    if isfield(opts,'minFunc'), opts = getopts(oldOpts,opts.minFunc); else, opts = oldOpts; end
    mode = opts.mode;
    nlZ_v = zeros(opts.repeat,1);
    dnlZ_v = zeros(opts.repeat,1);
    hyp_new = cell(opts.repeat,1);
    if mode < 4
        for i = 1:opts.repeat
            switch mode
                case 1
                    oldOpts.showIts = 1;
                    opts = getopts(oldOpts,opts);
                    if isempty(varargin)
                        [hyp_new{i}, nlZ_v(i), dnlZ_v(i)] = gradient_descent(hyp,@gp,opts,inf_,mean_,cov_,lik_,x,y);
                    else
                        [hyp_new{i}, nlZ_v(i), dnlZ_v(i)] = gradient_descent(hyp,@gp,opts,inf_,mean_,cov_,lik_,x,y, varargin{:});
                    end
                case 2
                    oldOpts.MaxFunEvals = 200;
                    oldOpts.Method = 'qnewton';
                    oldOpts.progTol = eps;
                    oldOpts.optTol = eps;
                    opts = getopts(oldOpts,opts);
                    [hyp_new{i}, nlZ_v(i), dnlZ_v(i)]=minimize_minfunc(hyp, @gp, opts, inf_, mean_, cov_, lik_, x, y);
                case 3
                    oldOpts.MaxFunEvals = 200;
                    opts = getopts(oldOpts,opts);
                    [hyp_new{i}, funval]=minimize(hyp, @gp, -opts.MaxFunEvals, inf_, mean_, cov_, lik_, x, y);
                    nlZ_v(i)=min(funval);
            end
            hyp = rand_hyp(opts.rndmInt);
        end
    else
        if ~isempty(varargin), algo_data = varargin{1}; end
        l= algo_data.l;
        direct.showits = 1;
        direct.maxevals = 1000;
        direct.maxits = 500;
        direct.maxdeep = 500;
        oldOpts.direct = direct;
        opts = getopts(oldOpts,opts);
        cond = cat(1,opts.rndmInt.mean,opts.rndmInt.lik,opts.rndmInt.cov([1;l],:));
        Problem.f=@(hyp_) DirectGP(hyp_, inf_, mean_, cov_, lik_, x, y, l,hyp);
        [~, hyp_] = Direct(Problem,cond,opts.direct);
        hyp_new=[];
        hyp_new.mean = hyp_(1);
        hyp_new.lik = log(hyp_(2));
        hyp_new.cov = hyp.cov;
        hyp_new.cov([end;l])=log(hyp_(3:end));
        [nlZ,dnlZ] = gp(hyp_new, inf_, mean_, cov_, lik_, x, y);
        return
    end
    [~,id] = min(nlZ_v);
    nlZ = nlZ_v(id);
    dnlZ = dnlZ_v(id);
    hyp_new = hyp_new{id};
    if nargout > 3
        i = id;
    end
end

function [hyp] = rand_hyp(rndmInt)
    hyp.cov = rndmInt.cov(:,1)+(rndmInt.cov(:,2)-rndmInt.cov(:,1)).*rand(size(rndmInt.cov,1),1);
    hyp.mean = rndmInt.mean(:,1)+(rndmInt.mean(:,2)-rndmInt.mean(:,1)).*rand(size(rndmInt.mean,1),1);
    hyp.lik = log(rndmInt.lik(:,1)+(rndmInt.lik(:,2)-rndmInt.lik(:,1)).*rand(size(rndmInt.lik,1),1));
end

function [nLz] = DirectGP(hyp__, inf_, mean_, cov_, lik_, x, y,l,hyp_)
    hyp.mean = hyp__(1);
    hyp.lik = log(hyp__(2));
    hyp.cov = hyp_.cov;
    hyp.cov([end;l]) = log(hyp__([3:end]));
    [nLz]=gp(hyp, inf_, mean_, cov_, lik_, x, y);
end
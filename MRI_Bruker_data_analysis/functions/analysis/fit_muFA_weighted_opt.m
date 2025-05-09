function out = fit_muFA_weighted_opt(opt)

out.lsqcurvefit = optimset('lsqcurvefit');

%increase the number of function evaluations for more accuracy
% out.lsqcurvefit.MaxFunEvals = 400;
% out.lsqcurvefit.MaxIter=400;
% out.lsqcurvefit.TolFun = 1e-006;
% out.lsqcurvefit.TolX = opt.fit.TolFun;
% out.lsqcurvefit.MaxFunEvals = 10000;
% out.lsqcurvefit.MaxIter=10000;
% out.lsqcurvefit.TolFun = 1e-008;
% out.lsqcurvefit.TolX = opt.fit.TolFun;
% out.lsqcurvefit = [];
out.lsqcurvefit.Display = 'off';

out.data_is_normalized = opt.data_is_normalized;
out.b_low_lim = opt.b_low_lim;

%filter with cut-off value b0 and steepness of filter s
out.weight.b0 = 2e9;
out.weight.s = 2;

% out.lb =   [0.5*s0_guess   1e-12  1e-26 1e-26];
% out.ub =   [1.5*s0_guess   3e-9   1e-12 1e-12];

out.lb =   [0.8   1e-12  1e-26 1e-26];
out.ub =   [1.2   3e-9   1e-12 1e-12];

if isfield(opt,'n_inits')
    out.n_inits = opt.n_inits;
end


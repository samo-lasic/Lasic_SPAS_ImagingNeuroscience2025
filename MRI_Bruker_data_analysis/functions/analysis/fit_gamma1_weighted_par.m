function Pout = fit_gamma1_weighted_par(b,sig,opt)

% use low b data to estimate s0
if ~opt.data_is_normalized
    %low_b_ind = find(abs(b-min(b)) < 1e6);
    low_b_ind = find(b < opt.b_low_lim);
    if isempty(low_b_ind)
        low_b_ind = [1 2]';
    end

    sig_low = sig(low_b_ind);
    b_low = b(low_b_ind);

    X = [ones(length(b_low),1) -b_low];
    P = X\log(sig_low);
    s0 = exp(P(1));
    %s0 = mean(sig(low_b_ind));

    opt.lb(1) = opt.lb(1) * s0;
    opt.ub(1) = opt.ub(1) * s0;
end
            

%------------------ weighting function for gamma fit ---------------
% filter with center 'c' and width 'w'
f = @(x,c,w) 1 - 1/2 * (tanh( (x-c) / w ) + 1);

% filter with cut-off value b0 and steepness of filter s
% f = @(x,b0,s) 1/(1+(b/b0).^(2*s));

weight = f(b/max(b),opt.weight.b0/max(b),opt.weight.s);

Pin =  (opt.lb + opt.ub)/2;
Pnorm = Pin;
Xnorm = mean(b);
Ynorm = mean(sig(:));
%[Pout,resnorm,residuals,flag] = lsqnonlin('fit_gamma1_error',...
Pout = Pnorm.*lsqnonlin('fit_gamma1_error',...
    Pin./Pnorm,opt.lb./Pnorm,opt.ub./Pnorm,opt.lsqcurvefit,b(:)/Xnorm,sig/Ynorm,weight,Pnorm,Xnorm,Ynorm);



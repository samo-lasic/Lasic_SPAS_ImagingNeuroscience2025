function Pout = fit_gamma1_par(b, sig, opt)

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

%------------------ bounds ---------------
Pin =  (opt.lb + opt.ub)/2;
Pnorm = Pin;
Xnorm = mean(b);
Ynorm = mean(sig(:));
%[Pout,resnorm,residuals,flag] = lsqnonlin('fit_gamma1_error',...
% Pout = lsqnonlin('fit_gamma1_error',...
%     Pin./Pnorm,opt.bounds.lb./Pnorm,opt.bounds.ub./Pnorm,opt.fit,b/Xnorm,sig/Ynorm,weight,Pnorm,Xnorm,Ynorm);

Pout = Pnorm.*lsqcurvefit('fit_gamma1',Pin./Pnorm,b(:)/Xnorm,sig/Ynorm,opt.lb./Pnorm,opt.ub./Pnorm,opt.lsqcurvefit,Pnorm,Xnorm,Ynorm);



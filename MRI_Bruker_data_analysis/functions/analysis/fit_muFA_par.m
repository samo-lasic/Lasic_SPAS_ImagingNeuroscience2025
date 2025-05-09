function Pout = fit_muFA_par(b,sig,opt)

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

    if ~isnan(s0)
        opt.lb(1) = opt.lb(1) * s0;
        opt.ub(1) = opt.ub(1) * s0;
    end
end

%---------------------------------
Xnorm = mean(b);
Ynorm = mean(sig(:));
%[Pout,resnorm,residuals,flag] = lsqnonlin('fit_gamma1_error',...
% Pout = lsqnonlin('fit_gamma1_error',...
%     Pin./Pnorm,opt.bounds.lb./Pnorm,opt.bounds.ub./Pnorm,opt.fit,b/Xnorm,sig/Ynorm,weight,Pnorm,Xnorm,Ynorm);

if isfield(opt,'n_inits')
    RSS = Inf;
    for n = 1:opt.n_inits

        Pin = opt.lb + rand(size(opt.lb)) .* (opt.ub - opt.lb); % guess
        Pnorm = Pin;

        Pout_tmp = Pnorm.*lsqcurvefit('fit_muFA',Pin./Pnorm,b(:)/Xnorm,sig/Ynorm,opt.lb./Pnorm,opt.ub./Pnorm,opt.lsqcurvefit,Pnorm,Xnorm,Ynorm);

        sig_fit = fit_muFA(Pout_tmp,b);
        RSS_tmp = sum((sig_fit-sig).^2); % residual sum of squares

        if RSS_tmp < RSS
            Pout = Pout_tmp;
            RSS = RSS_tmp;
        end
    end

else
    Pin =  (opt.lb + opt.ub)/2;
    Pnorm = Pin;

    Pout = Pnorm.*lsqcurvefit('fit_muFA',Pin./Pnorm,b(:)/Xnorm,sig/Ynorm,opt.lb./Pnorm,opt.ub./Pnorm,opt.lsqcurvefit,Pnorm,Xnorm,Ynorm);
end


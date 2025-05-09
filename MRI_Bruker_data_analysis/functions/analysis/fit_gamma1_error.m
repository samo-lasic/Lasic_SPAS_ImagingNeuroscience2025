%% function that reports the errorfit_gamma1fit_gamma1
function error = fit_gamma1_error(Pin,Xin,Yin,W,Pnorm,Xnorm,Ynorm);


% weight the error according to the |WEIGHT| vector

error = (fit_gamma1(Pin,Xin,Pnorm,Xnorm,Ynorm)-Yin).*W;

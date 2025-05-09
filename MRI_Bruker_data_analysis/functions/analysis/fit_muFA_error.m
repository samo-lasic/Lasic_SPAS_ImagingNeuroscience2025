%% function that reports the error
function error = fit_muFA_error(Pin,Xin,Yin,W,Pnorm,Xnorm,Ynorm);


% weight the error according to the |WEIGHT| vector

error = (fit_muFA(Pin,Xin,Pnorm,Xnorm,Ynorm)-Yin).*W;

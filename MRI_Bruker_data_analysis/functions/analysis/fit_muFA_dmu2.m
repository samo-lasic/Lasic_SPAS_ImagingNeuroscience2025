function Y = fit_muFA_dmu2(Pin,Xin,Pnorm,Xnorm,Ynorm);

if nargin == 2
    Pnorm = ones(size(Pin));
    Xnorm = 1;
    Ynorm = 1;
end

Pin = Pin.*Pnorm;
Xin = Xin*Xnorm;

Y0 = Pin(1);
D = Pin(2);
mu2_iso = Pin(3);
mu2_aniso = Pin(4);



npoints = length(Xin)/2;
Xin = Xin(1:npoints);

Y_aniso = Y0.*(1 + Xin.*mu2_aniso./D).^(-D.^2./mu2_aniso);
Y_iso = Y0.*(1 + Xin.*mu2_iso./D).^(-D.^2./mu2_iso);


Y = [Y_aniso; Y_iso]./Ynorm;

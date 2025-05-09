function Y = fit_gamma1(Pin,Xin,Pnorm,Xnorm,Ynorm)

if nargin == 2
    Pnorm = ones(size(Pin));
    Xnorm = 1;
    Ynorm = 1;
end

Pin = Pin.*Pnorm;
Xin = Xin*Xnorm;

Y0 = Pin(1);
D = Pin(2);
mu2 = Pin(3);

Y = Y0*((1 + Xin.*mu2/D).^(-1*D^2/mu2));
Y = Y/Ynorm;

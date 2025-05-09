function img = img_gamma_correction(img, gamma_cor)
A = max(abs(img(:)));
img = img ./ A;

%     p = 2; figure(2), clf, hold on
%     histogram(img(:)), histogram(img(:).^p)
%     p = 2; figure(2), clf, hold on
%     histogram(img(:)), histogram(sin(img(:)*pi/2).^p)
%     p = 2; figure(2), clf, hold on
%     histogram(img(:)), histogram(tanh(img(:)/p)/tanh(1/p))

if gamma_cor.type == 1
    img = img.^gamma_cor.p;
elseif gamma_cor.type == 2
    img = sin(img * pi/2).^gamma_cor.p;
elseif gamma_cor.type == 3
    img = tanh(img/gamma_cor.p)/tanh(1/gamma_cor.p);
end
img = img * A;
img = img * gamma_cor.A;


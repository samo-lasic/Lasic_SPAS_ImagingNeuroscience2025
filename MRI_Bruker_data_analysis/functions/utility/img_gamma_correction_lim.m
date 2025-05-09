function img = img_gamma_correction_lim(img, LIM, p)
% change contrast within the LIM range

ind = find(img(:) >= LIM(1) & img(:) <= LIM(2));
MIN = min(img(ind));
RANGE = max(img(ind)) - MIN;

% map to range 0-1
img(ind) = (img(ind) - MIN) / RANGE;

img(ind) = img(ind).^p;

% back to image range
img(ind) = MIN + RANGE * img(ind);




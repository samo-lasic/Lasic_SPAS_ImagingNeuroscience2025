function [g, b] = help_scale_gradient_from_b_range(bmin, bmax, Nb, b_at_100, p)
% returns gradient scales (in percent) for log-spaced b-values (if p == 0) else the values are power spaces
% absolute b-values don't matter, only the ratio bmin/bmax

if nargin < 5
    p = 0;
end

if p == 0

    b = logspace(log10(bmin),log10(bmax),Nb);

else
    b = linspace(0,1,Nb);
    b = bmin + (bmax-bmin) * b.^p;
end

g = round(100*sqrt(b/b_at_100));

b = round(b);



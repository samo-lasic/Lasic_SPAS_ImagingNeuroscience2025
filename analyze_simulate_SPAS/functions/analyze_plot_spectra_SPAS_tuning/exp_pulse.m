function val = exp_pulse(val, plato, par)

% Exponential shape defined by plateau, band, and p as parameters:
% - plateau = 0: __/\__ 
% - plateau = 1: --\/--


band = par.band;
p = par.p;
amp = par.amp;

if band > 0
    val = amp*exp(-1/band*val.^p);
else
    val = ones(size(val));
end

if plato == 1
    val = 1-val;
end
end

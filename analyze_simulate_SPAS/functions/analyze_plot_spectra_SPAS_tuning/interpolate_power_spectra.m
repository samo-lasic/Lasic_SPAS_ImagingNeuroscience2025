function [f_interp, PS_interp] = interpolate_power_spectra(f,PS,Ninterp)
Ninterp = 2*round(Ninterp/2)+1;
for i = 1:3
    for j = 1:3
        PS_interp(:,i,j) = interp1(1:length(PS),PS(:,i,j),linspace(1,length(PS),Ninterp));
    end
end

f_interp = interp1(1:length(f),f,linspace(1,length(f),Ninterp));
end
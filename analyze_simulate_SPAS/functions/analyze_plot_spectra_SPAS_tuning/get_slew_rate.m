function slew_rate = get_slew_rate(g,u,dt)
for n = 1:size(u,1)
    gu = sqrt(3)*g*u(n,:)';
    slew_rate.v(n) = max(abs(diff(gu)/dt));
end
% mean value, std, relative std
slew_rate.m = mean(slew_rate.v);
slew_rate.std = std(slew_rate.v);
slew_rate.rstd = slew_rate.std/slew_rate.m;
end



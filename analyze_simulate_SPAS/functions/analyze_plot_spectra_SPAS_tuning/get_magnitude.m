function mag = get_magnitude(g,u,dt)
for n = 1:size(u,1)
    gu = sqrt(3)*g*u(n,:)';
    mag.v(n) = max(abs(gu));
end
% mean value, std, relative std
mag.m = mean(mag.v);
mag.std = std(mag.v);
mag.rstd = mag.std/mag.m;
end

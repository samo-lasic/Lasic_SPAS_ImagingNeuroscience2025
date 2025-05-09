function C = color_map_RGB(x)
% x from 0 to 1
r = 1-2*x;
r(r<0) = 0;
g = 2*x;
g(g>1) = 2 - g(g>1);
b = 2*x - 1;
b(b<0) = 0;
%plot(x,b)

C = [r g b];
end
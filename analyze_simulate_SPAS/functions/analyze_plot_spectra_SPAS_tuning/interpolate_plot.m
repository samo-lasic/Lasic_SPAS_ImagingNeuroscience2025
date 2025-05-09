
function [X1, Y1] = interpolate_plot(X,Y,N)
if N > 0
    X1 = linspace(min(X),max(X),N);
    Y1 = interp1(X,Y,X1);
else
    X1 = X;
    Y1 = Y;
end
end
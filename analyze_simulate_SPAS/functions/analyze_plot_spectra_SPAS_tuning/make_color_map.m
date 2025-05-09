
function col = make_color_map(val,Ncolors)
% colormap for vector val
if Ncolors == 2 % 2 colors
    col(:,1) = (1-val);
    col(:,2) = val*0;
    col(:,3) = val;
elseif Ncolors == 1 % 3 colors
    tmp = 1-2*val; tmp(tmp<0) = 0;
    col(:,1) = tmp;
    tmp = 2*val; tmp(tmp>1) = 2-tmp(tmp>1);
    col(:,2) = tmp;
    tmp = -1+2*val; tmp(tmp<0) = 0;
    col(:,3) = tmp;
else
    col = ones(length(val),3);
end
end


function Tensor = tensor_info(DT)
% outputs [V, l, m, FA, Delta, Eta]

[V,L] = eig(DT);
m = trace(DT)/3;
l = [L(1,1) L(2,2) L(3,3)];
[l index] = sort(l,'descend');
V = V(:,index);
%maxl = max(l);

tmp = l/m-1;
[~, index] = sort(abs(tmp),'descend');
Delta = (tmp(index(1)) - (tmp(index(2))+tmp(index(3)))/2)/3;
tmp = tmp(index(2))-tmp(index(3));
if abs(tmp) < 1e-9
    Eta = 0;
else
    Eta = abs(tmp/2/Delta);
end
Delta = round(100*Delta)/100;
Eta = round(100*Eta)/100;

tmp = V;
for i = 1:3
    if (tmp(1,i)<0 & tmp(2,i)>=0 & tmp(3,i)<0) | ...
            (tmp(1,i)<0 & tmp(2,i)<0 & tmp(3,i)>=0) | ...
            (tmp(1,i)<0 & tmp(2,i)<0 & tmp(3,i)<0)  | ...
            (tmp(1,i)>=0 & tmp(2,i)<0 & tmp(3,i)<0)

        tmp(:,i) = -tmp(:,i);
    end
end
V = tmp;
FA = sqrt(3/2*sum((l-m).^2)/sum(l.^2));

Tensor.V = V;
Tensor.l = l;
Tensor.m = m;
Tensor.FA = FA;
Tensor.Delta = Delta;
Tensor.Eta = Eta;


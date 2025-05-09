
function Tensor = tensor_shape(Ts_in)
% size(Ts_in) = 3 x 3 or n x 3 x 3

if numel(size(Ts_in)) == 2
    Ts(1,:,:) = Ts_in;
else
    Ts = Ts_in;
end

for n = 1:size(Ts,1)
    T = squeeze(Ts(n,:,:));
    [V,L] = eig(T);
    m = trace(T)/3;
    l = [L(1,1) L(2,2) L(3,3)];
    [l, index] = sort(l,'descend');
    V = V(:,index);

    if m == 0
        tmp = l-1; % to fix the m = 0 singularity
    else
        tmp = l/m-1;
    end
    [~, index] = sort(abs(tmp),'descend');
    Delta = (tmp(index(1)) - (tmp(index(2))+tmp(index(3)))/2)/3;
    tmp = tmp(index(2))-tmp(index(3));
    if abs(tmp) < 1e-9
        Eta = 0;
    else
        Eta = abs(tmp/2/Delta);
    end

    Tensor(n).Delta = round(100*Delta)/100;
    Tensor(n).Eta = round(100*Eta)/100;
    Tensor(n).FA = sqrt(3/2*sum((l-m).^2)/sum(l.^2));

    Tensor(n).V = V;
    Tensor(n).l = l;
    Tensor(n).m = m;

end
end

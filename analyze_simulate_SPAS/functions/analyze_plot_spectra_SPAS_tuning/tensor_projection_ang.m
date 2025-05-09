function Tij = tensor_projection_ang(T,theta,phi)
% T(ind1,ind2,ind3)
% ind1 = time or frequency
% ind2,ind3 = T indices
% u(n,i) are n directions with i=1...3 components

u = [sin(theta)*cos(phi) sin(theta)*sin(phi) cos(theta)];

Tij(:) = 0;
for ind2 = 1:3
    for ind3 = 1:3
        Tij = Tij + squeeze(T(:,ind2,ind3))*u(ind2)*u(ind3);
    end
end


end
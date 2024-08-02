function [up,vp] = applyBC(ni,p)
    np = size(p,1);                        % Example change number 6
    vp = zeros(np,1);
    up = zeros(np,1);
    for i = 1:np
        vp(i) = nod2dof(ni,p(i,1),p(i,2));      % DOF index
        up(i) = p(i,3);                         % DOF value
    end
end
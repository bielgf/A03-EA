function f = pointLoads(data,f,F)
    for i = 1:size(F,1)
        I = nod2dof(data.ni,F(i,1),F(i,2));
        f(I) = F(i,3); % add point loads to the global forces vector
    end

end

% Example change number 5
function [xel,Sel,Mbel,Mtel] = internalforcesFunction(nel,ni,nne,x,Tn,Td,Kel,u)

    xel = zeros(2,nel);                     % Initialize vectors
    Sel = zeros(2,nel);
    Mbel = zeros(2,nel);
    Mtel = zeros(2,nel);
    
    for e = 1:nel
        xel(:,e) = x(Tn(e,:));                  % Element's vertices coordinates
        uel = zeros(nne*ni,1);
        for ii = 1:nne*ni
            uel(ii) = u(Td(e,ii));              % Assign DOF value
        end
        fel_int = Kel(:,:,e)*uel;               % Compute internal forces
        
        Sel(1,e) = -fel_int(1);                 % Assign cross-section loads
        Sel(2,e) = fel_int(4);
        Mbel(1,e) = -fel_int(2);
        Mbel(2,e) = fel_int(5);
        Mtel(1,e) = -fel_int(3);
        Mtel(2,e) = fel_int(6);
    end
end
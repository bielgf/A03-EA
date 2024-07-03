function [sigma,s] = normalStressDistribution(x_prim,Tn,x_0_prim,y_0_prim,I_xx_prim,I_yy_prim,I_xy_prim,M_x_prim,M_y_prim)

    n_el = size(Tn,1);
    
    s = zeros(2,n_el);
    sigma = zeros(2,n_el);
    
    %Equivalent inertias
    I_x = I_xx_prim-(I_xy_prim)^2/(I_yy_prim);
    I_y = I_yy_prim-(I_xy_prim)^2/(I_xx_prim);
    
    %Equivalent moments
    M_x = M_x_prim + (M_y_prim*I_xy_prim)/I_yy_prim;
    M_y = M_y_prim + (M_x_prim*I_xy_prim)/I_xx_prim;
    
    for i = 1:n_el
    
        %Element's length
        delta_x_prim = x_prim(Tn(i,2),1)-x_prim(Tn(i,1),1);
        delta_y_prim = x_prim(Tn(i,2),2)-x_prim(Tn(i,1),2);
        le = sqrt((delta_x_prim)^2+(delta_y_prim)^2);
    
        %Arc-length position corresponding to second node
        if i > 1
            s(1,i) = s(2,i-1);
        end
        s(2,i) = s(1,i) + le;
    
        %Stress corresponding to second node
        sigma(1,i) = (x_prim(Tn(i,1),2) - y_0_prim)*M_x/I_x - (x_prim(Tn(i,1),1) - x_0_prim)*M_y/I_y;
        sigma(2,i) = (x_prim(Tn(i,1),2) - y_0_prim)*M_x/I_x - (x_prim(Tn(i,1),1) - x_0_prim)*M_y/I_y;
    end

end
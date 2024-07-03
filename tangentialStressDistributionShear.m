% Example change number 3

function [tau_s,s] = tangentialStressDistributionShear (x_prim,Tn,m,Tm,x_0_prim,y_0_prim,I_xx_prim,I_yy_prim,I_xy_prim,S_x_prim,S_y_prim,x_s_prim,y_s_prim,A_in,open)

    n_el = size(Tn,1);
    
    s = zeros(2,n_el);
    tau_s = zeros(2,n_el);
    
    q = 0;
    
    if open == 0
        q = (S_y_prim*(x_0_prim - x_s_prim) - S_x_prim*(y_0_prim - y_s_prim))/(2*A_in);
    end
    
    %Equivalent inertias and shear
    I_x = I_xx_prim - (I_xy_prim)^2/I_yy_prim;
    I_y = I_yy_prim - (I_xy_prim)^2/I_xx_prim;
    S_x = S_x_prim - (S_y_prim*I_xy_prim)/I_xx_prim;
    S_y = S_y_prim - (S_x_prim*I_xy_prim)/I_yy_prim;
    
    for e = 1:n_el
        t = m(Tm(e),1);
    
        delta_x_prim = x_prim(Tn(e,2),1) - x_prim(Tn(e,1),1);
        delta_y_prim = x_prim(Tn(e,2),2) - x_prim(Tn(e,1),2);
        le = sqrt((delta_x_prim)^2 + (delta_y_prim)^2);
    
        if e > 1
            s(1,e) = s(2,e-1);
        end
        s(2,e) = s(1,e) + le;

        tau_s(1,e) = q/t;
        q = q - (S_x*t*le*(delta_x_prim/2 + x_prim(Tn(e,1),1) - x_0_prim))/I_y - (S_y*t*le*(delta_y_prim/2 + x_prim(Tn(e,1),2) - y_0_prim))/I_x;
        tau_s(2,e) = q/t;   
    end

end
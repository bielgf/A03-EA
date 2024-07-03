function [tau_t,s] = tangentialStressDistributionTorsion(x_prim,Tn,m,Tm,M_z_prim,J,open,A_in)

    n_el = size(Tn,1);
    
    s = zeros(2,n_el);
    tau_t = zeros(2,n_el);
    
    for ii = 1:n_el
        t = m(Tm(ii),1);
    
        delta_x_prim = x_prim(Tn(ii,2),1) - x_prim(Tn(ii,1),1);
        delta_y_prim = x_prim(Tn(ii,2),2) - x_prim(Tn(ii,1),2);
        le = sqrt((delta_x_prim)^2 + (delta_y_prim)^2);
    
        if ii > 1
            s(1,ii) = s(2,ii-1);
        end
        
        s(2,ii) = s(1,ii) + le;
    
        tau_t(1,ii) = M_z_prim*t/J;
        tau_t(2,ii) = M_z_prim*t/J;
    
        if open == 0
            tau_t(1,ii) = M_z_prim/(2*A_in*t);
            tau_t(2,ii) = M_z_prim/(2*A_in*t);
        end
    end
    
end
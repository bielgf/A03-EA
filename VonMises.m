function [max_vms,mcp_sect,mcp_beam] = VonMises(data,mD1,xnod,Sel,Mbel,Mtel)
    
    [x_prim,Tm,Tn] = GeometricDiscret(data);
    
    n_el = size(Tn,1);

    max_vms = 0;                                                      %Maximum Von Mieses stress across all elements
    mcp_sect = [1*10^(10),1*10^(10)];                                 %Most critical point of the section
    mcp_beam = 1*10^(10);                                             %Most critial point of the beam
%     max_vms_sect = zeros(data.nel,1);                                 %Maximum Von Mieses stress across the section
    
    for ii = 1:data.nel
        [x_0_prim,y_0_prim,x_s_prim,y_s_prim,~,I_xx_prim,I_yy_prim,I_xy_prim,J,A_in] = sectionProperties(data,x_prim,Tn,mD1,Tm,data.open);
        [sigma,~] = normalStressDistribution(x_prim,Tn,x_0_prim,y_0_prim,I_xx_prim,I_yy_prim,I_xy_prim,-Mbel(1,ii),0); 
        [tau_s,~] = tangentialStressDistributionShear(x_prim,Tn,mD1,Tm,x_0_prim,y_0_prim,I_xx_prim,I_yy_prim,I_xy_prim,0,Sel(1,ii),x_s_prim,y_s_prim,A_in,data.open); 
        [tau_t,~] = tangentialStressDistributionTorsion(x_prim,Tn,mD1,Tm,Mtel(1,ii),J,data.open,A_in); 
    
        for jj = 1:n_el                                               
            sigma_vm = sqrt(sigma(1,jj)^2+3*(tau_t(1,jj)+tau_s(1,jj))^2);
            
            if sigma_vm > max_vms
                max_vms = sigma_vm;
                mcp_sect = x_prim(jj,:);
                mcp_beam = xnod(ii);
            end
    
            % if sigma_vm > max_vms_sect(ii)
            %     max_vms_sect(ii) = sigma_vm;    
            % end
        end
    end

end
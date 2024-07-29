% Change example number 2

function [x_0_prim,y_0_prim,x_s_prim,y_s_prim,A_tot,I_xx_prim,I_yy_prim,I_xy_prim,J,A_in] = sectionProperties(x_prim,Tn,m,Tm,open)

    n_el = size(Tn,1);
    
    le = zeros(n_el,1);
    A = zeros(n_el,1);
    hat_x_prim = zeros(n_el,2);
    
    A_tot = 0;
    x_0_prim = 0;
    y_0_prim = 0;
    
    for ii = 1:n_el
        %Element's thickness
        t = m(Tm(ii),1);
        
        %Element's area and centroid position
        delta_x_prim = x_prim(Tn(ii,2),1) - x_prim(Tn(ii,1),1);
        delta_y_prim = x_prim(Tn(ii,2),2) - x_prim(Tn(ii,1),2);
        le(ii) = sqrt((delta_x_prim)^2+(delta_y_prim)^2);
        A(ii) = t*le(ii);
        hat_x_prim(ii,:) = (x_prim(Tn(ii,2),:) + x_prim(Tn(ii,1),:))/2;
    
        %Cummulative area and centroid position:
        A_tot = A_tot + A(ii);
        x_0_prim = x_0_prim + A(ii)*hat_x_prim(ii,1);
        y_0_prim = y_0_prim + A(ii)*hat_x_prim(ii,2);
    end
    
    %Centroid
    x_0_prim = x_0_prim/A_tot;
    y_0_prim = y_0_prim/A_tot;
    
    
    I_xx_prim = 0;
    I_yy_prim = 0;
    I_xy_prim = 0;
    J = 0;
    A_in = 0;
    
    
    for ee = 1:n_el
        %Element's thickness
        t = m(Tm(ee),1);
    
        %Inertia with respect to centroid
        delta_x_prim = x_prim(Tn(ee,2),1) - x_prim(Tn(ee,1),1);
        delta_y_prim = x_prim(Tn(ee,2),2) - x_prim(Tn(ee,1),2);
        I_xx_prim = I_xx_prim + A(ee)*(delta_y_prim)^2/12 + A(ee)*(hat_x_prim(ee,2) - y_0_prim)^2;
        I_yy_prim = I_yy_prim + A(ee)*(delta_x_prim)^2/12 + A(ee)*(hat_x_prim(ee,1) - x_0_prim)^2;
        I_xy_prim = I_xy_prim + A(ee)*(delta_y_prim)*(delta_x_prim)/12 + A(ee)*(hat_x_prim(ee,1) - x_0_prim)*(hat_x_prim(ee,2) - y_0_prim);
        J = J + le(ee)*t^3/3; 
    
        if open == 0
            A_mid = cross([x_prim(Tn(ee,1),1) - x_0_prim , x_prim(Tn(ee,1),2) - y_0_prim , 0],[delta_x_prim , delta_y_prim , 0]);
            A_in = A_in + norm(A_mid)/2;
            J = J + le(ee)/t;
        end
    end
    
    if open == 0
        J = 4*A_in^2/J;
    end
    
    q1 = 0;
    q2 = 0;
    x_s_prim = x_0_prim;
    y_s_prim = y_0_prim;
    
    for jj = 1:n_el   
        %Element's thickness
        t = m(Tm(jj),1);
    
        %Coefficients
        delta_x_prim = x_prim(Tn(jj,2),1) - x_prim(Tn(jj,1),1);
        delta_y_prim = x_prim(Tn(jj,2),2) - x_prim(Tn(jj,1),2);
    
        A1 = (I_yy_prim*delta_y_prim/2 - I_xy_prim*delta_x_prim/2)/(I_xx_prim*I_yy_prim - (I_xy_prim)^2);
        B1 = (I_yy_prim*(x_prim(Tn(jj,1),2) - y_0_prim) - I_xy_prim*(x_prim(Tn(jj,1),1) - x_0_prim))/(I_xx_prim*I_yy_prim - (I_xy_prim)^2);
        
        A2 = (I_xx_prim*delta_x_prim/2 - I_xy_prim*delta_y_prim/2)/(I_xx_prim*I_yy_prim - (I_xy_prim)^2);
        B2 = (I_xx_prim*(x_prim(Tn(jj,1),1) - x_0_prim) - I_xy_prim*(x_prim(Tn(jj,1),2) - y_0_prim))/(I_xx_prim*I_yy_prim - (I_xy_prim)^2);
    
        C = (x_prim(Tn(jj,1),1) - x_0_prim)*delta_y_prim - (x_prim(Tn(jj,1),2) - y_0_prim)*delta_x_prim;
    
        %Moment contribution
        x_s_prim = x_s_prim + C*(q1 - t*le(jj)*(A1/3+B1/2));
        y_s_prim = y_s_prim + C*(q2 + t*le(jj)*(A2/3+B2/2));
    
        %Shear flow update
        q1 = q1 - t*le(jj)*(A1+B1);
        q2 = q2 + t*le(jj)*(A2+B2);
    end

    if open == 0
        x_s_prim = x_0_prim;
        y_s_prim = y_0_prim;
    end
end
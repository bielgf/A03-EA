function StudyOfConvergence(data,mD2,x_s_prim,pD2,section)
    
    el_conv = [4 8 16 32 64 128 256 512];
    u_conv = zeros(3,size(el_conv,2));
    
    for nel = 1:size(el_conv,2)
        data.nel = el_conv(1,nel);

        [xnod,fe,me] = GetForceMomentElement(data,x_s_prim);
    
        TnD2 = [1:data.nel ; 2:data.nel+1]';
    
        TmD2 = ones(data.nel,1);
        
        TdD2 = connectDOF(data,TnD2);

        FD2 = [
                find(xnod == data.be)       1       -data.Me*data.g
                find(xnod == data.be)       3       -data.Me*data.g*(((data.d + data.xi_p) - x_s_prim) - data.ze)
                ];

        data.nnod = size(xnod,2);
        data.ndof = data.nnod*data.ni;

        [Kel] = stiffnessFunction(data,xnod',TnD2,mD2,TmD2);
        [fel] = forceFunction(data,xnod',TnD2,fe,me);
        [K,F] = assemblyFunction(data,TdD2,Kel,fel);
        [up,vp] = applyBC(data,pD2);
        [F] = pointLoads(data,F,FD2);
        [u,~] = solveSystem(data,K,F,up,vp);

        u_conv(:,nel) = u(end-2:end);
    end
    
    error_uy = zeros(1,size(u_conv,2)-1);
    error_rotx = zeros(1,size(u_conv,2)-1);
    error_rotz = zeros(1,size(u_conv,2)-1);

    for ii = 1:size(u_conv,2)-1
        error_uy(1,ii) = abs((u_conv(1,ii) - u_conv(1,end))/(u_conv(1,end)))*100;
        error_rotx(1,ii) = abs((u_conv(2,ii) - u_conv(2,end))/(u_conv(2,end)))*100;
        error_rotz(1,ii) = abs((u_conv(3,ii) - u_conv(3,end))/(u_conv(3,end)))*100;
    end
    
    figure(12)
    semilogx(el_conv(1:end-1),error_uy,el_conv(1:end-1),error_rotx,el_conv(1:end-1),error_rotz)
    legend('error in u_y','error in rot_x','error in rot_z')
    grid on
    xlabel('Beam''s number of elements' )
    ylabel('Error (%)')
    title(sprintf('Study of convergence for %s section',section))
end
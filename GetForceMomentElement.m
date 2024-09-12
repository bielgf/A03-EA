function [xnod,fe,me] = GetForceMomentElement(nel,aeroP,lambda,geomP)
    
    rhoinf = aeroP.rhoInf;
    vinf   = aeroP.vInf;
    cl     = aeroP.Cl;
    g      = aeroP.g;
    
    b    = geomP.wingspan;
    c    = geomP.chord;
    d    = geomP.beamWidth;
    xi_p = geomP.chiP;
    za   = geomP.aeroCenter;
    zm   = geomP.gravCenter;
    xsp  = geomP.xShearCenter;

    % Geometric validation
    if mod(nel,4) ~= 0
        fprintf('Invalid number of elements for beam''s geometric discretization\n')
        fe = -1;
        me = -1;
        return
    end
    
    % Distance computing
    zeta_s = (d + xi_p) - xsp;
    
    dLift = zeta_s - za;
    dWeight = zeta_s - zm;

    % Node computing
    xnod = 0:(b/nel):b;
    
    % Force computing
    x = sym('x','real');
    lx = 1/2*rhoinf*(vinf^2)*c*cl*sqrt(1 - (x/b)^2);
    lx = matlabFunction(lx);
    lx = lx(xnod);

    wx(1:size(xnod,2)) = lambda*g;

    
    % Vector's initialization
    fe = zeros(1,nel);
    me = zeros(1,nel);
    
    for ii = 1:nel
        Lift_e = (lx(ii) + lx(ii+1))/2;
        W_e = (wx(ii) + wx(ii+1))/2;

        fe(ii) = Lift_e - W_e;

        M_Lift_e = Lift_e*dLift;
        M_W_e = -W_e*dWeight;

        me(ii) = M_Lift_e + M_W_e;
    end
end
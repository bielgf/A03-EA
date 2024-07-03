function [xnod,fe,me] = GetForceMomentElement(data,xsp)
    
    % Geometric validation
    if mod(data.nel,4) ~= 0
        fprintf('Invalid number of elements for beam''s geometric discretization\n')
        fe = -1;
        me = -1;
        return
    end
    
    % Data
    b = data.b;
    rhoinf = data.rhoinf;
    vinf = data.vinf;
    c = data.c;
    Cl = data.Cl;
    lambda = data.lambda;
    g = data.g;
    
    % Distance computing
    zeta_s = (data.d + data.xi_p) - xsp;
    
    dLift = zeta_s - data.za;
    dWeight = zeta_s - data.zm;

    % Node computing
    xnod = 0:(b/data.nel):data.b;
    
    % Force computing
    x = sym('x','real');
    lx = 1/2*rhoinf*(vinf^2)*c*Cl*sqrt(1 - (x/b)^2);
    lx = matlabFunction(lx);
    lx = lx(xnod);

    wx(1:size(xnod,2)) = lambda*g;

    
    % Vector's initialization
    fe = zeros(1,data.nel);
    me = zeros(1,data.nel);
    
    for ii = 1:data.nel
        Lift_e = (lx(ii) + lx(ii+1))/2;
        W_e = (wx(ii) + wx(ii+1))/2;

        fe(ii) = Lift_e - W_e;

        M_Lift_e = Lift_e*dLift;
        M_W_e = -W_e*dWeight;

        me(ii) = M_Lift_e + M_W_e;
    end
end
function fel = forceFunction(nne,ni,nel,x,Tn,fe_bar,me_bar)
    fel = zeros(nne*ni,nel);
    for ei = 1:nel
        le = abs(x(Tn(ei,2),1) - x(Tn(ei,1),1));                        % Get the beam's lenght
    
        f_b = fe_bar(1,ei)*le*[1/2   le/12   0   1/2     -le/12  0]';      % Contribution of distributed shear load
    
        f_t = me_bar(1,ei)*le*[0  0  1/2   0   0   1/2]';                  % Contribution of distributed torsion
    
        f_ast = f_b + f_t;
        fel(:,ei) = f_ast;
    end
end
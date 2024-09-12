function Kel = stiffnessFunction(nne,ni,nel,x,Tn,m,Tm)
    Kel = zeros(nne*ni,nne*ni,nel);
    for ei = 1:nel
        le = abs(x(Tn(ei,2),1) - x(Tn(ei,1),1));        % Get the beam's lenght
    
        E = m(Tm(ei),1);                                % Get the beam's properties
        G = m(Tm(ei),2);
        I = m(Tm(ei),3);
        J = m(Tm(ei),4);
        
        b.E  = E;
        b.I  = I;
        b.le = le;
        BendMatrix = BendingStiffMatrixAssembly(b);
        BendMatrix.assembleMatrix();
        K_b = BendMatrix.Kb;
    
        K_t = (G*J)/(le)*[                              % Stiffness matrix for torsion
                            0   0   0   0   0   0
                            0   0   0   0   0   0
                            0   0   1   0   0   -1
                            0   0   0   0   0   0
                            0   0   0   0   0   0
                            0   0   -1  0   0   1
                            ];
    
        K_ast = K_b + K_t;
        Kel(:,:,ei) = K_ast;
    end
end
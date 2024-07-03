function Kel = stiffnessFunction(data,x,Tn,m,Tm)
    Kel = zeros(data.nne*data.ni,data.nne*data.ni,data.nel);
    for ei = 1:data.nel
        le = abs(x(Tn(ei,2),1) - x(Tn(ei,1),1));        % Get the beam's lenght
    
        E = m(Tm(ei),1);                                % Get the beam's properties
        G = m(Tm(ei),2);
        I = m(Tm(ei),3);
        J = m(Tm(ei),4);
        
        K_b = (E*I)/(le^3)*[                            % Stiffness matrix for bending
                            12      6*le    0   -12     6*le    0
                            6*le    4*le^2  0   -6*le   2*le^2  0
                            0       0       0   0       0       0           
                            -12     -6*le   0   12      -6*le   0
                            6*le    2*le^2  0   -6*le   4*le^2  0
                            0       0       0   0       0       0
                            ];
    
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
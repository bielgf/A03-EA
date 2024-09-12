function Kel = stiffnessFunction(nne,ni,nel,x,Tn,m,Tm)
    Kel = zeros(nne*ni,nne*ni,nel);
    for ei = 1:nel
        b.le = abs(x(Tn(ei,2),1) - x(Tn(ei,1),1));        % Get the beam's lenght
    
        b.E = m(Tm(ei),1);                                % Get the beam's properties
        b.G = m(Tm(ei),2);
        b.I = m(Tm(ei),3);
        b.J = m(Tm(ei),4);

        BendMatrix = BendingStiffMatrixAssembly(b);
        BendMatrix.assembleMatrix();
        Kb = BendMatrix.Kb;
    
        TorMatrix = TorsionStiffMatrixAssembly(b);
        TorMatrix.assembleMatrix();
        Kt = TorMatrix.Kt;
    
        K = Kb + Kt;
        Kel(:,:,ei) = K;
    end
end
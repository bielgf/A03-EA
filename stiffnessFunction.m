function Kel = stiffnessFunction(nne,ni,nel,x,nodalConnec,beamProp,materialConnec)
    Kel = zeros(nne*ni,nne*ni,nel);
    for ei = 1:nel
        b.le = abs(x(nodalConnec(ei,2),1) - x(nodalConnec(ei,1),1));        % Get the beam's lenght
    
        b.E = beamProp(materialConnec(ei),1);                                % Get the beam's properties
        b.G = beamProp(materialConnec(ei),2);
        b.I = beamProp(materialConnec(ei),3);
        b.J = beamProp(materialConnec(ei),4);

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
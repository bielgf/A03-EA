classdef StiffnessFunctionClass < handle
    
    properties (Access = private)
        beamParams
        connec
        beamProp
    end

    properties (Access = public)
        Kel
    end

    methods (Access = public)
        
        function obj = StiffnessFunctionClass(cParams)
            obj.init(cParams);
        end

        function computeStiffnessMatrix(obj)
            obj.compute();
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.beamParams = cParams.beamP;
            obj.connec     = cParams.con;
            obj.beamProp   = cParams.beamPr;
        end

        function compute(obj)
            nne = obj.beamParams.numNodesElem;
            ni  = obj.beamParams.numDOFsNode;
            nel = obj.beamParams.numElements;
            x   = obj.beamParams.xGlobal';
            nC  = obj.connec.nodalConnec;
            mC  = obj.connec.materialConnec;
            bP  = obj.beamProp;

            obj.Kel = zeros(nne*ni,nne*ni,nel);
            for ei = 1:nel
                b.le = abs(x(nC(ei,2),1) - x(nC(ei,1),1));
            
                b.E = bP(mC(ei),1);
                b.G = bP(mC(ei),2);
                b.I = bP(mC(ei),3);
                b.J = bP(mC(ei),4);
        
                BendMatrix = BendingStiffMatrixAssembly(b);
                BendMatrix.assembleMatrix();
                Kb = BendMatrix.Kb;
            
                TorMatrix = TorsionStiffMatrixAssembly(b);
                TorMatrix.assembleMatrix();
                Kt = TorMatrix.Kt;
            
                K = Kb + Kt;
                obj.Kel(:,:,ei) = K;
            end
        end

    end

end
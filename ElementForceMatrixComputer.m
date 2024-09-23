classdef ElementForceMatrixComputer < handle
    
    properties (Access = private)
        beamParams
        connec
    end

    properties (Access = public)
        fel
    end

    methods (Access = public)

        function obj = ElementForceMatrixComputer(cParams)
            obj.init(cParams);
        end

        function computeForceVector(obj)
            obj.compute();
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.beamParams = cParams.beamP;
            obj.connec     = cParams.con;
        end

        function compute(obj)
            nne = obj.beamParams.numNodesElem;
            ni  = obj.beamParams.numDOFsNode;
            nel = obj.beamParams.numElements;
            x   = obj.beamParams.xGlobal';
            fe  = obj.beamParams.forceElem;
            me  = obj.beamParams.momentElem;
            Tn  = obj.connec.nodalConnec;

            fElem = zeros(nne*ni,nel);
            for ei = 1:nel
                le = abs(x(Tn(ei,2),1) - x(Tn(ei,1),1));
            
                fb = fe(1,ei)*le*[1/2   le/12   0   1/2     -le/12  0]';
            
                ft = me(1,ei)*le*[0  0  1/2   0   0   1/2]';
            
                f = fb + ft;
                fElem(:,ei) = f;
            end
            obj.fel = fElem;
        end

    end

end
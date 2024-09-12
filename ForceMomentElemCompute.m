classdef ForceMomentElemCompute < handle
    
    properties (Access = private)
        numDOFsNode
        numElements
        geoP
        aeroParams
        lambda
    end

    methods (Access = public)
        
        function obj = ForceMomentElemCompute(cParams)
            obj.init(cParams);
        end

        function [xGlobal,forceElem,momentElem,totalDOFs] = compute(obj)
            nElem = obj.numElements;
            geomP = obj.geoP;
            aeroP = obj.aeroParams;
            l     = obj.lambda;
            [xGlobal,forceElem,momentElem] = GetForceMomentElement(nElem,aeroP,l,geomP);
            nnod = size(xGlobal,2);
            totalDOFs = nnod*obj.numDOFsNode;
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.numDOFsNode  = cParams.numDOFsNode;
            obj.numElements  = cParams.numElements;
            obj.geoP         = cParams.geoP;
            obj.aeroParams   = cParams.aeroParams;
            obj.lambda       = cParams.lambda;
        end

    end

end
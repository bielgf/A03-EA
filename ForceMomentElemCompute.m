classdef ForceMomentElemCompute < handle
    
    properties (Access = private)
        ndof
        nelem
        geoParams
        aeroParams
        lambda
    end

    methods (Access = public)
        
        function obj = ForceMomentElemCompute(cParams)
            obj.init(cParams);
        end

        function [xGlobal,forceElem,momentElem,totalDOFs] = compute(obj)
            nElem = obj.nelem;
            geomP = obj.geoParams;
            aeroP = obj.aeroParams;
            l     = obj.lambda;
            [xGlobal,forceElem,momentElem] = GetForceMomentElement(nElem,aeroP,l,geomP);
            nnod      = size(xGlobal,2);
            totalDOFs = nnod*obj.ndof;
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.ndof       = cParams.numDOFsNode;
            obj.nelem      = cParams.numElements;
            obj.geoParams  = cParams.geoP;
            obj.aeroParams = cParams.aeroParams;
            obj.lambda     = cParams.lambda;
        end

    end

end
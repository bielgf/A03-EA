classdef ForceMomentElemCompute < handle
    
    properties (Access = private)
        geoParams
        aeroParams
        lambda
        ndof
        nelem
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
            obj.geoParams  = cParams.geoParams;
            obj.aeroParams = cParams.aeroParams;
            obj.lambda     = cParams.lambda;
            obj.ndof       = cParams.numDOFsNode;
            obj.nelem      = cParams.numElements;
        end

    end

end
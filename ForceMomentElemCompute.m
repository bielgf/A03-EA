classdef ForceMomentElemCompute < handle
    
    properties (Access = private)
        numDOFsNode
        numElements
        wingspan
        rhoInf
        vInf
        chord
        Cl
        lambda
        g
        beamWidth
        chiP
        aeroCenter
        centerOfMass
        x_s_prim
    end

    methods (Access = public)
        
        function obj = ForceMomentElemCompute(cParams)
            obj.init(cParams);
        end

        function [xnod,fe,me,ndof] = compute(obj)
            [xnod,fe,me] = GetForceMomentElement(obj.numElements,obj.wingspan,obj.rhoInf,obj.vInf,obj.chord,obj.Cl,obj.lambda,obj.g,obj.beamWidth,obj.chiP,obj.aeroCenter,obj.centerOfMass,obj.x_s_prim);
            nnod = size(xnod,2);
            ndof = nnod*obj.numDOFsNode;
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.numDOFsNode       = cParams.numDOFsNode;
            obj.numElements      = cParams.numElements;
            obj.wingspan        = cParams.wingspan;
            obj.rhoInf   = cParams.rhoInf;
            obj.vInf     = cParams.vInf;
            obj.chord        = cParams.chord;
            obj.Cl       = cParams.Cl;
            obj.lambda   = cParams.lambda;
            obj.g        = cParams.g;
            obj.beamWidth        = cParams.beamWidth;
            obj.chiP     = cParams.chiP;
            obj.aeroCenter       = cParams.aeroCenter;
            obj.centerOfMass       = cParams.centerOfMass;
            obj.x_s_prim = cParams.x_s_prim;
        end

    end

end
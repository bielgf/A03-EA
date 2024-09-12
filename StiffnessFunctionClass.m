classdef StiffnessFunctionClass < handle
    
    properties (Access = private)
        numNodesElem
        numDOFsNode
        numElements
        xGlobal
        nodalConnec
        beamProp
        materialConnec
    end

    methods (Access = public)
        
        function obj = StiffnessFunctionClass(cParams)
            obj.init(cParams);
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.numNodesElem   = cParams.numNodesElem;
            obj.numDOFsNode    = cParams.numDOFsNode;
            obj.numElements    = cParams.numElements;
            obj.xGlobal        = cParams.xGlobal';
            obj.nodalConnec    = cParams.nodalConnec;
            obj.beamProp       = cParams.beamProp;
            obj.materialConnec = cParams.materialConnec;
        end

    end

end
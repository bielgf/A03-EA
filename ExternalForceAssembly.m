classdef ExternalForceAssembly < handle
    
    properties (Access = private)
        geomParams
        aeroParams
        beamParams
    end

    properties (Access = public)
        externalForce
    end

    methods (Access = public)
        
        function obj = ExternalForceAssembly(cParams)
            obj.init(cParams);
        end

        function assembly(obj)
            obj.assemblyMatrix();
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.geomParams = cParams.geomParams;
            obj.aeroParams = cParams.aeroParams;
            obj.beamParams = cParams.beamParams;
        end

        function assemblyMatrix(obj)
            [F11,F21,F13,F23] = obj.MatrixComponents();
            obj.externalForce = [F11, 1, F13 ; F21, 3, F23];
        end

        function [F11,F21,F13,F23] = MatrixComponents(obj)
            xG   = obj.beamParams.xGlobal;
            xE   = obj.beamParams.xEngine;
            zE   = obj.beamParams.zEngine;
            eM   = obj.beamParams.engineMass;
            g    = obj.aeroParams.g;
            bw   = obj.geomParams.beamWidth;
            chiP = obj.geomParams.chiP;
            xSc  = obj.geomParams.xShearCenter;
            
            F11 = find(xG == xE);
            F21 = F11;
            F13 = -eM*g;
            F23 = -eM*g*(((bw + chiP) - xSc) - zE);
        end

    end

end
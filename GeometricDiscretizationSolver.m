classdef GeometricDiscretizationSolver < handle

    properties (Access = private)
        geomParams
    end

    methods (Access = public)
        
        function obj = GeometricDiscretizationSolver(cParams)
            obj.init(cParams);
        end

        function [xSection,materialSecConnec,nodalSecConnec] = compute(obj)
            N1        = obj.geomParams.N1;
            N2        = obj.geomParams.N2;
            N3        = obj.geomParams.N3;
            h1        = obj.geomParams.h1;
            h2        = obj.geomParams.h2;
            beamWidth = obj.geomParams.beamWidth;
            [xSection,materialSecConnec,nodalSecConnec] = GeometricDiscret(N1,N2,N3,beamWidth,h1,h2);
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.geomParams = cParams;
        end

    end

end
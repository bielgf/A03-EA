classdef GeometricDiscretizationSolver < handle

    properties (Access = private)
        N1
        N2
        N3
        beamWidth
        h1
        h2
    end

    methods (Access = public)
        
        function obj = GeometricDiscretizationSolver(cParams)
            obj.init(cParams);
        end

        function [x_prim,Tm,Tn] = compute(obj)
            [x_prim,Tm,Tn] = GeometricDiscret(obj.N1,obj.N2,obj.N3,obj.beamWidth,obj.h1,obj.h2);
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.N1 = cParams.N1;
            obj.N2 = cParams.N2;
            obj.N3 = cParams.N3;
            obj.beamWidth  = cParams.beamWidth;
            obj.h1 = cParams.h1;
            obj.h2 = cParams.h2;
        end

    end

end
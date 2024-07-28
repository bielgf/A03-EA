classdef MAIN < handle
    
    properties (Access = private)
        data
    end

    properties (Access = public)
        x_prim
        Tm
        Tn
    end

    methods (Access = public)
        
        function obj = MAIN(cParams)
            obj.init(cParams);
        end

        function computeGeoDiscret(obj)
            obj.geometricdiscretization();
        end

    end

    methods (Access = private)
        
        function init(obj,cParams)
            obj.data = cParams;
        end

        function geometricdiscretization(obj)
            [obj.x_prim,obj.Tm,obj.Tn] = GeometricDiscret(obj.data);
        end

    end

end
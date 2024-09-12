classdef TorsionStiffMatrixAssembly < handle
    
    properties (Access = private)
        shearM
        shearI
        l
    end

    properties (Access = public)
        Kt
    end

    methods (Access = public)
        
        function obj = TorsionStiffMatrixAssembly(cParams)
            obj.init(cParams);
        end

        function assembleMatrix(obj)
            obj.assemble;
        end

    end

    methods (Access = private)
        
        function init(obj,cParams)
            obj.shearM  = cParams.G;
            obj.shearI  = cParams.J;
            obj.l = cParams.le;
        end

        function assemble(obj)
            G  = obj.shearM;
            J  = obj.shearI;
            le = obj.l;
            obj.Kt = (G*J)/(le)*[
                            0   0   0   0   0   0
                            0   0   0   0   0   0
                            0   0   1   0   0   -1
                            0   0   0   0   0   0
                            0   0   0   0   0   0
                            0   0   -1  0   0   1
                            ];
        end

    end

end
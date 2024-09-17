classdef BendingStiffMatrixAssembly < handle
    
    properties (Access = private)
        youngM
        xInertia
        l
    end

    properties (Access = public)
        Kb
    end

    methods (Access = public)
        
        function obj = BendingStiffMatrixAssembly(cParams)
            obj.init(cParams);
        end

        function assembleMatrix(obj)
            obj.assemble;
        end

    end

    methods (Access = private)
        
        function init(obj,cParams)
            obj.youngM    = cParams.E;
            obj.xInertia  = cParams.I;
            obj.l         = cParams.le;
        end

        function assemble(obj)
            E  = obj.youngM;
            I  = obj.xInertia;
            le = obj.l;
            obj.Kb = (E*I)/(le^3)*[                            
                            12      6*le    0   -12     6*le    0
                            6*le    4*le^2  0   -6*le   2*le^2  0
                            0       0       0   0       0       0           
                            -12     -6*le   0   12      -6*le   0
                            6*le    2*le^2  0   -6*le   4*le^2  0
                            0       0       0   0       0       0
                            ];
        end

    end

end
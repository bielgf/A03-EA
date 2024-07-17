classdef IterativeSolver < Solver
    
    properties (Access = private)
        LHS
        RHS
    end

    properties (Access = public)
        x
    end

    methods (Access = public)
        
        function obj = IterativeSolver(LHS,RHS)
            obj.init(LHS,RHS);
        end

        function compute(obj)
            obj.SolveSystem();
        end

    end

    methods (Access = private)
        
        function init(obj,LHS,RHS)
            obj.LHS = LHS;
            obj.RHS = RHS;
        end

        function SolveSystem(obj)
            obj.x = pcg(obj.LHS,obj.RHS,1e-8,10000);
        end
    
    end

end
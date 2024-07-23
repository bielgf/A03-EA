classdef DirectSolver < handle
    
    properties (Access = private)
        LHS
        RHS
    end

    properties (Access = public)
        x
    end

    methods (Access = public)
        
        function obj = DirectSolver(LHS,RHS)
            obj.init(LHS,RHS);
        end

        function compute(obj)
            obj.solveSystem();
        end

    end

    methods (Access = private)
        
        function init(obj,LHS,RHS)
            obj.LHS = LHS;
            obj.RHS = RHS;
        end

        function solveSystem(obj)
            obj.x = obj.LHS\obj.RHS;
        end
    
    end

end
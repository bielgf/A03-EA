classdef DirectSolver < Solver
    
    properties
        x
    end

    methods (Access = public)
        
        function obj = DirectSolveSystem(obj,LHS,RHS)
            obj.x = LHS\RHS;
        end

    end
end
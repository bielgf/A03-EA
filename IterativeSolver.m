classdef IterativeSolver < Solver
    
    properties
        x
    end
    
    methods (Access = public)
        
        function obj = IterativeSolveSystem(obj,LHS,RHS)
            obj.x = pcg(LHS,RHS,1e-10);
        end

    end
end
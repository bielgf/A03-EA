classdef Solver
    
    properties
        solution
    end

    methods (Access = public)
        
        function obj = SolveSystem(obj,LHS,RHS,method)
            
            switch method
                case 'Direct'
                    DS = DirectSolver();
                    obj.solution = DS.DirectSolveSystem(LHS,RHS);
                case 'Iterative'
                    IS = IterativeSolver();
                    obj.solution = IS.IterativeSolveSystem(LHS,RHS);
            end
        
        end

    end
end
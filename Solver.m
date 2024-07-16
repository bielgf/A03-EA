classdef Solver < handle

    methods (Access = public, Static)
        
        function obj = create(obj,LHS,RHS,method)
            
            switch method
                case 'Direct'
                    obj = DirectSolver();
                    %obj.solution = DS.DirectSolveSystem(LHS,RHS);
                case 'Iterative'
                    obj = IterativeSolver();
                    %obj.solution = IS.IterativeSolveSystem(LHS,RHS);
            end
        
        end

    end
end
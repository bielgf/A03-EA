classdef Solver < handle

    methods (Access = public, Static)
        
        function obj = create(LHS,RHS,method)
            
            switch method
                case 'Direct'
                    obj = DirectSolver(LHS,RHS);
                case 'Iterative'
                    obj = IterativeSolver(LHS,RHS);
            end
        
        end

    end
end
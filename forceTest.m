classdef forceTest < matlab.unittest.TestCase
    
    properties (Access = private)
        F
        FGood
    end

    methods (Access = public)
        function obj = forceTest(F)
            obj.init(F);
        end
    end

    methods (Access = private)
        
        function init(obj,F)
            load("resultsGood.mat",'FGood')
            obj.F = F;
            obj.FGood = FGood;
        end

    end

    methods (Test)

        function compareSolution(testCase)
            error = norm(testCase.F(:) - testCase.FGood(:))/norm(testCase.FGood(:));
            testCase.verifyLessThanOrEqual(error,1e-10);
        end

    end

end
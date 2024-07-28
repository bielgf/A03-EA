classdef displacementsTest < matlab.unittest.TestCase
    
    properties (Access = private)
        u
        uGood
    end

    methods (Access = public)
        function obj = displacementsTest(u)
            obj.init(u);
        end
    end

    methods (Access = private)
        
        function init(obj,u)
            load("resultsGood.mat",'uGood')
            obj.u = u;
            obj.uGood = uGood;
        end

    end

    methods (Test)

        function compareSolution(testCase)
            error = norm(testCase.u(:) - testCase.uGood(:))/norm(testCase.uGood(:));
            testCase.verifyLessThanOrEqual(error,1e-10);
        end

    end

end
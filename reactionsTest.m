classdef reactionsTest < matlab.unittest.TestCase
    
    properties (Access = private)
        r
        rGood
    end

    methods (Access = public)
        function obj = reactionsTest(r)
            obj.init(r);
        end
    end

    methods (Access = private)
        
        function init(obj,r)
            load("resultsGood.mat",'rGood')
            obj.r = r;
            obj.rGood = rGood;
        end

    end

    methods (Test)

        function compareSolution(testCase)
            error = norm(testCase.r(:) - testCase.rGood(:))/norm(testCase.rGood(:));
            testCase.verifyLessThanOrEqual(error,1e-10);
        end

    end

end
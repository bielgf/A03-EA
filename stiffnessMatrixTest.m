classdef stiffnessMatrixTest < matlab.unittest.TestCase
    
    properties (Access = private)
        K
        KGood
    end

    methods (Access = public)
        function obj = stiffnessMatrixTest(K)
            obj.init(K);
        end
    end

    methods (Access = private)
        
        function init(obj,K)
            load("resultsGood.mat",'KGood')
            obj.K = K;
            obj.KGood = KGood;
        end

    end

    methods (Test)

        function compareSolution(testCase)
            error = norm(testCase.K(:) - testCase.KGood(:));
            testCase.verifyLessThanOrEqual(error,1e-10);
        end

    end

end
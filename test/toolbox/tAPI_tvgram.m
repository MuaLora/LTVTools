classdef tAPI_tvgram < matlab.unittest.TestCase
%% TVGRAM 
% 
% API test for tvgram - tests basic functionality
    
    methods(Test)
        function tvgram_basic(testCase)
            
            % <Comment> tvgram for computing controllability gramian
            % </Comment>
            Time = linspace(0,5,10)';
            AData = -5+0.1*Time.^2;
            A = tvmat(AData,Time);
            B = 1; C = 1; D=0;
            G = tvss(A,B,C,D);
            Wc = tvgram(G,'c');
            testCase.verifyNotEmpty(Wc,'Wc should not be empty');
            testCase.verifyClass(Wc,'tvmat','Wc should be of type ''tvmat''');
        
            % <Comment> ODE Options and ODE Solver can be specified using
            % tvodeOptions </Comment>
            tvgramopt = tvodeOptions;
            tvgramopt.OdeOptions = odeset('RelTol',1e-4,'AbsTol',1e-6);
            Wo = tvgram(G,'o',tvgramopt);
            testCase.verifyNotEmpty(Wo,'Wo should not be empty');
            testCase.verifyClass(Wo,'tvmat','Wo should be of type ''tvmat''');
        end
    end
end
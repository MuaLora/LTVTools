classdef tNUM_tvgram < matlab.unittest.TestCase
    %% TVGRAM
    %
    % Numerical test for tvgram - tests for numerical accuracy and
    % reproducibility
    
    methods(TestMethodTeardown)
        function methodCleanup(testCase) %#ok<MANU>
            % Clean up the workspace by closing all the figures
            close('all');
        end
    end
    
    methods(Test)       
        
        function tvgram_SISO_numerical(testCase)
            % <Comment> Use tvgram to compute controllability and
            % observability gramian for stable SISO system </Comment>
            Time = linspace(0,5,100)';
            AData = -5 + 2*sin(Time); % Negative values (Periodic Signal)
            A = tvmat(AData,Time);
            B = 2; C = 1; D=0;
            G = tvss(A,B,C,D);
            Wc_stable = tvgram(G,'c');
            Wo_stable = tvgram(G,'O');
            figure;plot(Wc_stable);title('Controllability Gramian (Stable System)');
            figure;plot(Wo_stable);title('Observability Gramian (Stable System)');
            
            % <Comment> Use tvgram to compute controllability and
            % observability gramian for unstable SISO system </Comment>
            AData = 5 + 2*sin(Time); % Positive values (Periodic Signal)
            A = tvmat(AData,Time);
            G = tvss(A,B,C,D);
            Wc_unstable = tvgram(G,'c');
            Wo_unstable = tvgram(G,'O');
            figure;plot(Wc_unstable);title('Controllability Gramian (Unstable System)');
            figure;plot(Wo_unstable);title('Observability Gramian (Unstable System)');
            
            % Numerical Match with baseline data
            load('bl_tvgram_siso.mat','Wc_stable_b','Wo_stable_b',...
                'Wc_unstable_b','Wo_unstable_b');
            testCase.verifyEqual(Wc_stable,Wc_stable_b,'AbsTol',1e-8,...
                'Gramian Values Mismatch with baseline data');
            testCase.verifyEqual(Wo_stable,Wo_stable_b,'AbsTol',1e-8,...
                'Gramian Values Mismatch with baseline data');
            testCase.verifyEqual(Wc_unstable,Wc_unstable_b,'AbsTol',1e-8,...
                'Gramian Values Mismatch with baseline data');
            testCase.verifyEqual(Wo_unstable,Wo_unstable_b,'AbsTol',1e-8,...
                'Gramian Values Mismatch with baseline data');
        end
    end
end
classdef tNUM_tvlqr < matlab.unittest.TestCase
    % Tests the TVLQR function for numerical accuracy
    
    methods(TestMethodTeardown)
        function methodCleanup(testCase) %#ok<MANU>
            % Clean up the workspace by closing all the figures
            close('all');
        end
    end
    
    methods(Test)
        function testLTI(testCase)
            % Kirk, D. E. (2012). Optimal control theory:
            % an introduction. Courier Corporation.(Ex 5.2-2) Second Order
            % LTI system
            
            % System Parameters
            A = [0 1;2 -1];
            B = [0; 1];
            C = eye(2);
            D = 0;
            Q = [2 0;0 1];
            R = 1/2;
            sys_ss = ss(A,B,C,D);
            sys = tvss(sys_ss);
            F = zeros(2);
            
            % Use TVLQR to integrate CRDE
            [K,P] = tvlqr(sys,Q,R,F,[0 15]);
            
            % Reshape the data
            Pdata1 = reshape(P.Data(1,1,:),size(P.Time));
            Pdata12 = reshape(P.Data(1,2,:),size(P.Time));
            Pdata21 = reshape(P.Data(2,1,:),size(P.Time));
            Pdata2 = reshape(P.Data(2,2,:),size(P.Time));
            Pdata = [Pdata1 Pdata12 Pdata21 Pdata2];
            
            Kdata1 = reshape(K.Data(1,1,:),size(P.Time));
            Kdata12 = reshape(K.Data(1,2,:),size(P.Time));
            Kdata = [Kdata1 Kdata12];
            
            % Plot the data
            figure;
            plot(P.Time,Pdata);
            legend('P11','P12','P21','P22');
            title('Solution of Riccati Equation')
            % NOTE: The P plot should match with Page no. 218, Figure 5-9
            % (a)
            figure;plot(K.Time,Kdata);
            legend('K1','K12');
            title('State Feedback Gain Values');
            
            % Compare with inf horizon lqr with AbsTol of 1e-4
            [K_IH,S_IH] = lqr(sys_ss,Q,R); 
            testCase.verifyEqual(K_IH,Kdata(1,:),'AbsTol',1e-4,...
                'K data do not match with inf horizon data');
            testCase.verifyEqual(S_IH,reshape(Pdata(1,:),size(S_IH)),'AbsTol',1e-4,...
                'Riccati Sol. P data do not match with inf horizon data');
            
            % Compare with baseline
            load('bl_tvlqr_lti.mat','PbData','KbData');
            
            % Verify that results do match with baseline
            testCase.verifyEqual(Kdata,KbData,'AbsTol',1e-8,...
                'K data do not match with baseline');
            testCase.verifyEqual(Pdata,PbData,'AbsTol',1e-8,...
                'P data do not match with baseline');
            
            % Do simulation and plot state trajectories
            x0 = [-4 4];
            Acl = A - B*K;
            Bcl = sys.B*0;
            R = tvmat(zeros(length(Acl.Time),1),Acl.Time); % Track 0
            sys_cl = tvss(Acl,Bcl,sys.C,sys.D);
            X = tvlsim(sys_cl,R,x0);
            Xdata1 = reshape(X.Data(1,:,:),[length(X.Time) 1]);
            Xdata2 = reshape(X.Data(2,:,:),[length(X.Time) 1]);
            Xdata = [Xdata1 Xdata2];
            
            % Control inputs
            Xeval = evalt(X,K.Time,'tvmat');
            U = -K*Xeval;
            Udata = reshape(U.Data,size(K.Time));
            
            % Plot Page no. 218, Figure 5-9 (b)
            figure;
            plot(X.Time,Xdata); hold on;
            plot(U.Time,Udata);
            legend('X1*','X2*','U*');
            title('Optimal Control Input and State Trajectories');
            
            % Compare with baseline
            load('bl_tvlqr_lti.mat','XbData','UbData');
            
            % Numerical Verification
            testCase.verifyEqual(Xdata,XbData,'AbsTol',1e-8,...
                'X data do not match with baseline');
            testCase.verifyEqual(Udata,UbData,'AbsTol',1e-8,...
                'U data do not match with baseline');
        end
    end
end
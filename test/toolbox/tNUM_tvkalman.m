classdef tNUM_tvkalman < matlab.unittest.TestCase
    % Tests the TVKALMAN function for numerical accuracy
    
    methods(TestMethodTeardown)
        function methodCleanup(testCase) %#ok<MANU>
            % Clean up the workspace by closing all the figures
            close('all');
        end
    end
    
    methods(Test)
        
        function testSISO(testCase)
            % Specify Time Grid
            Tf = 30;
            Time = linspace(0,Tf,1000)';
            
            % System Under Test Parameters and time varying system
            %      .
            %      x = A(t)x(t) + B(t)u(t) + w(t)       {State equation}
            %      y(t) = C(t)x(t) + D(t)u(t) + v(t)    {Measurements}
            A = -3;
            AData = A*ones(size(Time));
            Amat = tvmat(AData,Time);
            B = 1;
            C = 1;
            D = 0;
            G = 1;
            H = 0;
            sys_tvss_kf = tvss(Amat,[B G],C,[D H]);
            sys_tvss = tvss(Amat,B,C,D);
            
            % For Inf. Horizon analysis
            %      .
            %      x = Ax + Bu + Gw            {State equation}
            %      y = Cx + Du + Hw + v        {Measurements}            
            sys_ss = ss(A,[B G],C,[D H]);
            
            % Initial Conditions
            P0 = 1;
            x0 = -1;
            
            % Process Noise Covariance
            Qn = 0.01;
            
            % Measurement Noise Covariance
            Rn = 0.1;
            
            % Make sure LTV system is Observable
            Wo = tvgram(sys_tvss_kf,'o');
            figure;plot(Wo);title('Observability Gramian (Stable System)');
            
            % Use TVKALMAN to integrate CRDE
            [Kest,L,P] = tvkalman(sys_tvss_kf,Qn,Rn,P0,[0 Tf]);
            %[Ae,Be,Ce,De] = ssdata(Kest); %#ok<ASGLU>
            
            % Plot the data
            figure;plot(P);
            title('Solution of Riccati Equation (Error Covariance P)');
            figure;plot(L);
            title('Kalman Filter Gain Values');
            
            % Use Infinite Horizon kalman filter to compute gains and
            % covariances (Ricatti Diffrential Equation solutions)
            [Kest_IH,L_IH,P_IH] = kalman(sys_ss,Qn,Rn); %#ok<ASGLU>
            testCase.verifyEqual(L_IH,L.Data(end),'AbsTol',1e-4,...
                'Kalman Gain L data do not match with inf horizon data');
            testCase.verifyEqual(P_IH,P.Data(end),'AbsTol',1e-4,...
                'Riccati Sol. P data do not match with inf horizon data');
            
            % ========================
            % Create Baseline MAT file
            % ========================
            LData = reshape(L.Data,size(L.Time));
            PData = reshape(P.Data,size(P.Time));
            % LbData = LData;PbData = PData;Wo_b=Wo;
            % save('bl_tvkalman_siso.mat','LbData','PbData','Wo_b');
            
            % Compare with baseline
            load('bl_tvkalman_siso.mat','LbData','PbData','Wo_b');
            
            % Verify that results do match with baseline
            testCase.verifyEqual(LData,LbData,'AbsTol',1e-8,...
                'L data do not match with baseline');
            testCase.verifyEqual(PData,PbData,'AbsTol',1e-8,...
                'P data do not match with baseline');
            testCase.verifyEqual(Wo,Wo_b,'AbsTol',1e-8,...
                'Wo data do not match with baseline');
            
            % Evaluate Kalman Filter Performance
            [A,B,C,D] = ssdata(sys_tvss);
            Aaug = [A      zeros(size(A-L*C));...
                    L*C    A-L*C];
            w_zeros = tvmat(zeros(size(sys_tvss_kf.Time)),sys_tvss_kf.Time);
            w_ones = tvmat(ones(size(sys_tvss_kf.Time)),sys_tvss_kf.Time);
            v_zeros = tvmat(zeros(size(sys_tvss_kf.Time)),sys_tvss_kf.Time);
            v_ones = tvmat(ones(size(sys_tvss_kf.Time)),sys_tvss_kf.Time);
            Baug = [B w_ones  v_zeros;...
                B w_zeros v_zeros];
            Caug = [C C*0];
            Daug = [D w_zeros v_ones];
            aug_sys = tvss(Aaug,Baug,Caug,Daug);
            
            % Create Input Signals
            t = sys_tvss_kf.Time;
            n = length(sys_tvss_kf.Time);
            W = tvmat(sqrt(Qn)*randn(n,1),t);
            V = tvmat(sqrt(Rn)*randn(n,1),t);
            % Control input is 0, initial condition response
            u = tvmat(zeros(n,1),t);
            U = [u;W;V];
            
            % Augmented Plant initial condition
            X0 = [x0;x0+4];
            
            % Simulate System
            [~,X] = tvlsim(aug_sys,U,X0);
            X = evalt(X,t,'tvmat');
            
            % Plot Results
            figure;plot(X);title('Kalman Filter State Estimate Results');
            legend('x','xhat');
            
            % ========================
            % Append Baseline MAT file
            % ========================
            % Xb = X; % Augmented State Vector
            % save('bl_tvkalman_siso.mat','Xb','-append');
            
            % Compare with Baseline
            load('bl_tvkalman_siso.mat','Xb');
            testCase.verifyEqual(X.Data,Xb.Data,'AbsTol',1e-1,'RelTol',1e-1,...
                'X data do not match with baseline');
            testCase.verifyEqual(X.Data(:,:,end),[0;0],'AbsTol',1e-1,...
                'RelTol',1e-1,'State Eventually does not go to zero');
        end
    end
end
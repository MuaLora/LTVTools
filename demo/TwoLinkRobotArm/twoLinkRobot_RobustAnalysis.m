%% twoLinkRobot_RobustAnalysis
% Analysis for Linearized Two-Link Robot Model
%
% This file analyzes the robust performance of a two-link robot arm along 
% a  pre-specified trajectory. A summary of the robot and corresponding 
% dynamics can be found in Section 2.3 of Ref [1]. The modeling and 
% analysis code in this file has been updated from original code written 
% by Moore [2].
%
% Ref:
% [1] R. Murray, Z. Li, and S. Sastry. A Mathematical Introduction to 
% Robot Manipulation. CRC Press, 1994.
% http://www.cds.caltech.edu/~murray/books/MLS/pdf/mls94-manipdyn_v1_2.pdf
%
% [2] R. Moore. Finite horizon robustness analysis using integral quadratic
% constraints. Master’s thesis, University of California, Berkeley, 2015.

%% Construct linearized model and compute state feedback
% See twoLinkRobot_BuildLTVModel for details on the nonlinear model,
% linearized model, and state-feedback law.
twoLinkRobot_BuildLTVModel;
close(1:3);

%% Parameters for Robustness Analysis 

% IQC Parameters: Pole p<0 and Filter order v>=0
v = 1; 
p = -10;  

% LMI Time Grid
Nt = 20; 
tlmi = linspace(0,Tf,Nt);

% Spline Basis Time Grid
Nsp = 10; 
tSp = linspace(0,Tf,Nsp);

% Uncertainty Norm Bound
DelNorm = 0.8;

%% Robust Open-Loop Analysis
if true
    % Uncertain OL System:  Gunc = Fu(Gnom,Delta)
    % Delta is unit norm-bounded LTI uncertainty.  Gnom is constructed 
    % so that uncertainty enters into second input channel as 0.8*Delta.
    Gnom = [0 0 sqrt(DelNorm); G(1:2,:)*[0 1 0; sqrt(DelNorm) 0 1] ];
    
    % Robust (worst-case) L2 to E Gain
    % Expected 995 (MS Thesis)
    fprintf('\n----------------- Open-Loop Analysis\n');
    [gOL,wcinfoOL] = tvrobL2toE(Gnom,v,p,Tf,tlmi,tSp);
end

%% Robust Closed-Loop Analysis

% Uncertain CL System:  Tunc = Fu(Tnom,Delta)
systemnames = 'G K';
inputvar = '[w; d1; d2]';
outputvar = '[d2-K(2) ; G(1:2)]';
input_to_G = '[d1-K(1); d2-K(2)+w]';
input_to_K = '[G]';
cleanupsysic = 'yes';
Tnom = sysic;

Lscl = diag([sqrt(DelNorm) 1 1]);
Rscl = diag([sqrt(DelNorm) 1 1]);
Tnom = Lscl*Tnom*Rscl;

% Robust (worst-case) L2 to E Gain
% Expected 0.062 (MS Thesis) 
fprintf('\n----------------- Closed-Loop Analysis\n');
[gCL,wcinfoCL] = tvrobL2toE(Tnom,v,p,Tf,tlmi,tSp);

%% Evaluate closed-loop with a specific bad perturbation

% Wrap in "Bad" Perturbation
% The *.Mat file contains 100 random samples of perturbations (with
% norm ||Delta||=0.8). The code below selects a perturbation and 
% normalizes it to have ||Delta||=1.
% [Note: The selected Delta is the one that maximizes the closed-loop gain
%  out of these 100 random samples.  The calculation of the closed-loop
%  norms on these 100 samples is skipped to reduce computation.]
load RandomDelta.mat;
DeltaBad = DeltaSim{11} / norm(DeltaSim{11},inf);
Tbad  = lft( DeltaBad ,Tnom);

% Evaluate closed-loop gain of Tbad
% Gain lower bound should be close to 0.0575
fprintf('\n----------------- Evaluate "Bad" Perturbation \n');
[gWC,dWC] = tvnorm(Tbad,[0 Tf],'L2toE');
fprintf('\n Closed-loop gain with worst-case Delta = %4.4f',gWC(2))

% Simulate linear system and evaluate gain
yBad = tvlsim(Tbad,dWC);
gWC2 = norm( evalt(yBad,Tf) ) / tvnorm(dWC);
fprintf('\n Closed-loop gain with worst-case Delta and dist. = %4.4f\n',gWC2)

%% Simulate and Plot Trajectories / Norm Bound

% Disturbances are scaled to have norm ||d|| = dL2norm
dL2norm = 5;

% e(T) ball: Center at final (trim) angles with radius gCL*dL2norm
rBall = gCL*dL2norm;
cBall = evalt( etabar(1:2), Tf);
NBall = 50;
aBall = linspace(0,2*pi,NBall);
th1Ball = cBall(1) + rBall*cos(aBall) + 2*pi;
th2Ball = cBall(2) + rBall*sin(aBall);

% Convert e(T) ball to Cartesian coordinates
x2Ball = L1*cos(th1Ball) + L2*cos( th1Ball + th2Ball );
y2Ball = L1*sin(th1Ball) + L2*sin( th1Ball + th2Ball );

% Draw e(T) ball and corresponding Cartesian coordinate region
figure(5)
clf;
axis([0 5, -3.5 0.5]);
patch(th1Ball,th2Ball,'c');
hold on

figure(6)
clf
axiswidth = L1 + L2 + 0.2*(L1 + L2);
axis equal; 
axis(axiswidth*[-1, 1, -1 ,1]);hold on;
patch(x2Ball,y2Ball,'c');
hold on;

% Load Random Disturbances and append worst-case disturbance
load RandomDisturbances.mat;
d = [d, {dWC}];

% Simulate linearized system 
for i = 1:numel(d)
    % Display count i
    if floor(i/10)==ceil(i/10)
        fprintf('\n i=%d ',i)
    end
    
    % Disturbances are scaled to have norm ||d||=dL2norm    
    di = d{i};
    di = di*dL2norm/tvnorm(di);

    % Simulate linear system
    yi = tvlsim(Tbad,di);

    % Superimpose the trim trajectory to obtain actual angles
    % (See twoLinkRobot_BuildLTVModel for trim trajectory construction)    
    theta = evalt( etabar(1:2), yi.Time, 'tvmat');
    theta.InterpolationMethod = 'Linear';
    theta = yi+ theta;
    theta1 = theta(1)+2*pi;
    theta2 = theta(2);

    % Plot trajectory and mark final point
    figure(5)
    plot(theta1,theta2);    
    theta1f = evalt(theta1,Tf);
    theta2f = evalt(theta2,Tf);
    plot(theta1f,theta2f,'ko','MarkerFaceColor','w');
        
    % Cartesian coordinates of link 2 tip
    x2 = L1*cos(theta1) + L2*cos(theta1+theta2); 
    y2 = L1*sin(theta1) + L2*sin(theta1+theta2);
    
    % Plot trajectory and mark final point
    figure(6)
    plot(x2,y2);
    x2f = evalt(x2,Tf);
    y2f = evalt(y2,Tf);
    plot(x2f,y2f,'ko','MarkerFaceColor','w');        
end
fprintf('\n')

% Add labels and legends
figure(5)
grid on
xlabel('\theta_1 (rads)'); 
ylabel('\theta_2 (rads)'); 
title(['Closed-Loop with ||d|| <= ' num2str(dL2norm)]); 
hold off

figure(6)
grid on
xlabel('x (m)'); ylabel('y (m)'); 
title(['Closed-Loop with ||d|| <= ' num2str(dL2norm)]); 
hold off



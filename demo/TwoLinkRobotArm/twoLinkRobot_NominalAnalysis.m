%% twoLinkRobot_NominalAnalysis
% Analysis for Linearized Two-Link Robot Model
%
% This file analyzes the nominal performance of a two-link robot arm along 
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

%% Nominal Open-Loop Analysis

% Nominal Open Loop L2 to E Gain
% Expected 528.9 (MS Thesis) 
[gOLNom,dOL,infoOL] = tvnorm(G(1:2,:),[T0 Tf],'L2toE');
fprintf('\n Bounds on Nominal OL Gain = [%4.4f, %4.4f] \n',...
    gOLNom(2),gOLNom(1));
dOL = 0.001*dOL/tvnorm(dOL);

%% Open-Loop Linear and Nonlinear Simulations

% Linear Simulation with Worst-Case Disturbance
[yLin,xLin]=tvlsim(G(1:2,:),dOL);
yf = evalt(yLin,yLin.Time(end));
gLin = norm(yf) / tvnorm(dOL);
fprintf(' Nominal OL Gain Linear Sim = %4.4f\n',gLin);
etaLin = etabar;
etaLin.InterpolationMethod = 'Linear';
etaLin = etaLin +  evalt(xLin,Tgrid,'tvmat');

% Nonlinear Simulation
d = dOL;
sim('TwoLinkRobotOL',[T0 Tf]);
etaNL = tvmat(eta',tsim);
etaNL = evalt(etaNL,Tgrid,'tvmat');

% Plot Results
figure(1); clf
subplot(311)
plot(dOL.Time,dOL(1),'b',dOL.Time,dOL(2),'r--');
ylabel('d (Nm)');
legend('d1','d2');
title('Simulations with ||dOL||=0.001')

subplot(312)
plot(Tgrid,etabar(1),'k-.', Tgrid,etaNL(1),'r', ...
    Tgrid,etaLin(1),'b--');
ylabel('\theta_1 (rad)');
legend('Trim','NL','Lin');

subplot(313)
plot(Tgrid,etabar(2),'k-.', Tgrid,etaNL(2),'r', ...
    Tgrid,etaLin(2),'b--');
xlabel('Time (sec)');
ylabel('\theta_2 (rad)');

figure(2); clf
plot(etabar(1),etabar(2),'k-.', etaNL(1),etaNL(2),'r', ...
    etaLin(1),etaLin(2),'b--');
xlabel('\theta_1 (rad)');
ylabel('\theta_2 (rad)');
legend('Trim','NL','Lin');
title('Simulations with ||dOL||=0.005')


%% Nominal Closed-Loop Analysis

% Build closed-loop from disturbance to linearized state x
T = feedback(G,K);

% Nominal Closed-Loop Loop L2 to E Gain
% Expeced 0.05509 (MS Thesis)
[gCLNom,dCL,infoCL] = tvnorm(T(1:2,:),[T0 Tf],'L2toE');
fprintf('\n Bounds on Nominal CL Gain = [%4.4f, %4.4f] \n',...
    gCLNom(2),gCLNom(1));
dCL = 20*dCL/tvnorm(dCL);

%% Closed-Loop Linear and Nonlinear Simulations

% Linear Simulation with Worst-Case Disturbance
[yLin,xLin]=tvlsim(T(1:2,:),dCL);
yf = evalt(yLin,yLin.Time(end));
gLin = norm(yf) / tvnorm(dCL);
fprintf(' Nominal CL Gain Linear Sim = %4.4f\n',gLin);
etaLin = etabar;
etaLin.InterpolationMethod = 'Linear';
etaLin = etaLin +  evalt(xLin,Tgrid,'tvmat');

% Nonlinear Simulation
d = dCL;
sim('TwoLinkRobotCL',[T0 Tf]);
etaNL = tvmat(eta',tsim);
etaNL = evalt(etaNL,Tgrid,'tvmat');

% Plot Results
figure(3); clf
subplot(311)
plot(dCL.Time,dCL(1),'b',dCL.Time,dCL(2),'r--');
ylabel('d (Nm)');
legend('d1','d2');
title('Simulations with ||dCL||=20')

subplot(312)
plot(Tgrid,etabar(1),'k-.', Tgrid,etaNL(1),'r', ...
    Tgrid,etaLin(1),'b--');
ylabel('\theta_1 (rad)');
legend('Trim','NL','Lin');

subplot(313)
plot(Tgrid,etabar(2),'k-.', Tgrid,etaNL(2),'r', ...
    Tgrid,etaLin(2),'b--');
xlabel('Time (sec)');
ylabel('\theta_2 (rad)');

figure(4); clf
plot(etabar(1),etabar(2),'k-.', etaNL(1),etaNL(2),'r', ...
    etaLin(1),etaLin(2),'b--');
xlabel('\theta_1 (rad)');
ylabel('\theta_2 (rad)');
legend('Trim','NL','Lin');
title('Simulations with ||dCL||=20')

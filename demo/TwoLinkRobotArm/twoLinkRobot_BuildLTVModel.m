%% twoLinkRobot_BuildLTVModel
% Build LTV model for Two-Link Robot 
%
% This file builds an LTV model for a two-link robot arm along a  pre-
% specified trajectory. A state-feedback law is also computed using 
% finite-horizon LQR. A summary of the robot and corresponding dynamics 
% can be found in Section 2.3 of Ref [1]. The modeling and analysis code 
% in this file has been updated from original code written by Moore [2].
%
% Ref:
% [1] R. Murray, Z. Li, and S. Sastry. A Mathematical Introduction to 
% Robot Manipulation. CRC Press, 1994.
% http://www.cds.caltech.edu/~murray/books/MLS/pdf/mls94-manipdyn_v1_2.pdf
%
% [2] R. Moore. Finite horizon robustness analysis using integral quadratic
% constraints. Master’s thesis, University of California, Berkeley, 2015.


%% Nonlinear Model parameters
% The configuration of the two-link robot is described by the joint
% angles [theta1; theta1].  Motors provide torques, tau=[tau1; tau2], 
% at the base of each link.  The nonlinear model has the form.
%     d(eta)/dt = f(eta,tau)
% where eta = [theta1;theta2; theta1dot; theta2dot]. Further details are
% provided in references [1] and [2].

% Robot model parameters
m1 = 3;         % Mass of link 1, kg
m2 = 2;         % Mass of link 2, kg
L1 = 0.3;       % Length of link 1, meters
L2 = 0.3;       % Length of link 2, meters
r1 = L1/2;      % Length to centroid of link 1, meters
r2 = L2/2;      % Length to centroid of link 2, meters
Iz1 = m1*L1^2/3; % Moment of inertia of link 1, kg*m^2
Iz2 = m2*L2^2/3; % Moment of inertia of link 2, kg*m^2

beta = m2*L1*r2;
delta = Iz2 + m2*r2^2;
alpha = Iz1 + Iz2 + m1*r1^2 + m2*(L1^2 + r2^2);

% Store key robot parameters in a structure
RobotParam.beta = beta;
RobotParam.delta = delta;
RobotParam.alpha = alpha;
RobotParam.L1 = L1;
RobotParam.L2 = L2;

%% Equilbrium State Trajectory and Input

% Time Horizon
T0 = 0;
Tf = 5;

% Specify (trim) trajectory for end of link 2 in Cartesian coordinates 
% Trajectory is specified as a cubic spline with a small number of interp.
% points.  The spline is then evaluated on a dense time grid.
x2Data = [-0.5 -0.5 0 0.5 0.5];
x2Time = Tf*[0 0.001 0.5 0.999 1];
x2bar = tvmat(x2Data,x2Time,'Spline');

Tgrid = linspace(T0,Tf,400);
x2bar = evalt(x2bar,Tgrid,'tvmat');
y2bar = 0.25*sin(5*x2bar)+0.1;

% Solve inverse kinematics problem for joint angles, joint velocities, 
% and torque inputs corresponding to the motion of link 2 end point.
% Use joint angles to compute cartesian points for end of link 1.
% The state is eta=[theta1; theta2; theta1_dot; theta2_dot] and etabar
% is the trim trajectory (in link angle coordinates).
[etabar,taubar] = twoLinkInverseKinematics(x2bar,y2bar,RobotParam);
x1bar = L1*cos(etabar(1));
y1bar = L1*sin(etabar(1));

% Plot equilibrium input torques
figure(1); clf
plot(Tgrid,taubar(1),'b',Tgrid,taubar(2),'r--','LineWidth',2);
legend('tau1','tau2')
ylabel('Torque (Nm)')
xlabel('Time (s)')
ylim([-12,10])
grid on;

% Simulate nonlinear system with equilibrium inputs
d = tvmat([0;0],[T0;Tf]);
eta0 = evalt(etabar,T0);
etaf = evalt(etabar,Tf);
sim('TwoLinkRobotOL',[T0 Tf]);
etasim = tvmat(eta',tsim);

x1sim = L1*cos(etasim(1)); 
y1sim = L1*sin(etasim(1));
x2sim = x1sim + L2*cos(etasim(1)+etasim(2)); 
y2sim = y1sim + L2*sin(etasim(1)+etasim(2));

% Plot Trajectory: Link Angles
% Note: Simulated results should match analytical trim (bar) results
figure(2); clf
plot( etabar(1), etabar(2), 'b', etasim(1), etasim(2), 'r--', ...
    eta0(1), eta0(2),'bo',etaf(1), etaf(2),'bx',...
    'LineWidth',3,'MarkerSize',10)
xlabel('\theta_1 (rad)');
ylabel('\theta_2 (rad)');
grid on;
set(gca, 'LooseInset', get(gca, 'TightInset'));

%%
% Plot Trajectory: Cartesian Coordinates and Snapshots of Robot Position
figure(3); clf
plot( x2bar, y2bar, 'b', x2sim, y2sim, 'r--', ...
    evalt(x2bar,T0), evalt(y2bar,T0),'bo', ...
    evalt(x2bar,Tf), evalt(y2bar,Tf),'bx', ...
    'LineWidth',3,'MarkerSize',10)
xlabel('x (m)');
ylabel('y (m)');
grid on;

hold on;
for ti = [1.5 2.5 4 5.0]    
    xi = evalt([x1bar;x2bar;y1bar; y2bar],ti);
    plot([0; xi(1:2)],[0; xi(3:4)],'k','LineWidth',3);
    plot(xi(1),xi(3),'ko','MarkerSize',8,'MarkerFaceColor','b');
    plot(xi(2),xi(4),'ko','MarkerSize',8,'MarkerFaceColor','r');
    plot(0,0,'ko','MarkerSize',8,'MarkerFaceColor','k');
end
hold off;
axis equal;
set(gca, 'LooseInset', get(gca, 'TightInset'));

text(0.01,-0.3,'t=1.8','FontSize',12);
text(-0.39,0.08,'t=2.5','FontSize',12);
text(0,0.3,'t=4.0','FontSize',12);
text(0.2,0.18,'t=5.0','FontSize',12);


%% Build LTV System
% States of G are x=eta-etabar.
% These states are perturbations around the trim trajectory etabar.
io(1) = linio('TwoLinkRobotOL/Input Torque',1,'input');
io(2) = linio('TwoLinkRobotOL/Two Link Robot Arm',1,'output');
sys = linearize('TwoLinkRobotOL',io,Tgrid);
G = tvss(sys,Tgrid);
[A,B,C,D] = ssdata(G);

%% Time-Varying LQR State-Feedback
Q = 100*diag([1 1 0.1 0.1]);
R = 0.1*diag([1,1]);
F = diag([1 1 0.1 0.1]);
K = tvlqr(A,B,Q,R,[],F,[T0,Tf]);

% Evaluate K on same grid as plant model G
K = evalt(K,Tgrid,'tvmat');


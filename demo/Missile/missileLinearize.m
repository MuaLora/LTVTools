%% missileLinearize

%% Model Parameters
% States are x:=[alpha; q; V]
% Input is delta_P
% Outputs are y = [x; nZ]
rho = 0.0023769;        % slugs/ft^3
d = 0.5;                % ft
A = 0.25*pi*d^2;        % ft^2
m = 250;                % slugs
Iyy = 49.4;             % slugs*ft^2


%% Simulate Nonlinear System and Plot Results
% The results will be treated as the trim trajectory.

% Specify IC and input for simulation
Tf = 4;
x0bar = [0.098;0;3050];
tbar = linspace(0,Tf,100);
tbar = tbar(:);
deltabar = -1+tbar/Tf;
deltabar = deltabar(:);

% Run simulation
delta_in = [tbar deltabar];
x0 = x0bar;
sim('firstSim',tbar);
ybar = ysim;

% Plot Results
figure(1); clf
plot(tsim,delta_sim);
ylabel('delta');
xlabel('Time (sec)')

figure(2); clf
for i=1:4
    subplot(2,2,i)
    plot(tsim,ybar(:,i))
    grid on;
end
subplot(2,2,1); ylabel('alpha');
subplot(2,2,2); ylabel('q');
subplot(2,2,3); ylabel('V'); xlabel('Time (sec)');
subplot(2,2,4); ylabel('nZ'); xlabel('Time (sec)');

%% Linearize along a trajectory
% Specify I/O points using LinIO.  The syntax is
%    IO=linio('blockname',PORTNUM,TYPE)
% The 'blockname" and "PORTNUM" specify a signal for either the 
% linearization input or output. I obtained the name by clicking on 
% the block and then typing "gcb" at the Matlab prompt.
% The "TYPE" input is either the linearization "input" or "output".
%
% Specify grid of "snapshot" times. The linearization is performed by
% simulating the system with the given inputs/ICs and linearizing along
% the trajectory at the specified snapshot times.  If the snapshot times
% are omitted then I think it simply performs the linearization at t=0.
Tgrid = linspace(0,Tf,25);
io(1) = linio('firstSim/From Workspace',1,'input');
io(2) = linio('firstSim/S-Function',1,'output');
G = linearize('firstSim',io,Tgrid);
[Ag,Bg,Cg,Dg]=ssdata(G);

%% Compare nonlinear and linearized responses

% Specify small change to input and IC
deltaPert = 0.2*sin(tbar);
deltaPert = deltaPert(:);
x0Pert = [0;0;10];

% Run nonlinear simulation
delta_in = [tbar deltabar+deltaPert];
x0 = x0bar+x0Pert;
sim('firstSim',tbar);

% Run linearized simulation
% I have most of the infrastructure already coded to do this in LTVtools....
fh = @(t,x) Glin(t,x,deltaPert,tbar,Ag,Bg,Tgrid);
[tlin,xlin] = ode45(fh,tbar,x0Pert);

% Plot Results
figure(3); clf
for i=1:4
    subplot(2,2,i)
    plot(tbar,ybar(:,i),'k'); hold on;
    plot(tbar,ysim(:,i),'b','LineWidth',2); hold on;
    if i<4
        plot(tbar,ybar(:,i)+xlin(:,i),'r--','LineWidth',2)
    end
    grid on;
end
subplot(2,2,1); ylabel('alpha');
legend('Trim','NL','Lin')
subplot(2,2,2); ylabel('q');
subplot(2,2,3); ylabel('V'); xlabel('Time (sec)');
subplot(2,2,4); ylabel('nZ'); xlabel('Time (sec)');

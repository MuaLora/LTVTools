function [a1,a2,a3,a4] = missileSfun(t,x,u,message,x0,rho,d,A,m,Iyy)

if message==0
   tmp = simsizes;
   tmp.NumContStates = 3; %[alpha;q;V]
   tmp.NumDiscStates = 0;
   tmp.NumOutputs = 4;  % [alpha,q,V,nz]
   tmp.NumInputs = 1; % [deltaP]
   tmp.DirFeedthrough = 1; % does y depend on u?
   tmp.NumSampleTimes = 1;
   a1 = simsizes(tmp);
   a2 = x0;  % x(0), initial condition, nx1 column
   a3 = [];  % always,
   a4 = [0 0]; % Sample time, offset
elseif message==1
   alpha = x(1);
   q = x(2);
   V = x(3);
   deltaP = u;
   CNz = -21*alpha - 24.5*alpha^3 + (-2+0.1*sqrt(abs(alpha))*6/pi)*deltaP;
   Cdap = 0.5*(alpha+deltaP)^2;
   CD = 0.2;
   CMy = -15.75*alpha - 36.75*alpha^3 +(-2+0.1*sqrt(abs(alpha))*42/pi)*deltaP;
   CMtd = -720;
   Fx = -0.5*rho*V^2*A*(Cdap + CD);
   Fz = 0.5*rho*V^2*A*CNz;
   My = 0.5*rho*V^2*A*d*(CMy + CMtd*d/2/V*q);
   ca = cos(alpha);
   sa = sin(alpha);
   ta = tan(alpha);
   Qab = 1/ca;
   alphaDot = ca^2*(Fz*Qab/m/V) + q - sa*ca*Fx*Qab/m/V;
   qDot = My/Iyy;
   VDot = 1/m/Qab*(Fx + Fz*ta);
   a1 = [alphaDot;qDot;VDot];
elseif message==3
   % outputs are [States;nZ]
   alpha = x(1);
   V = x(3);
   deltaP = u;
   CNz = -21*alpha - 24.5*alpha^3 + (-2+0.1*sqrt(abs(alpha))*6/pi)*deltaP;
   Fz = 0.5*rho*V^2*A*CNz;
   a1 = [x;Fz/m];
end

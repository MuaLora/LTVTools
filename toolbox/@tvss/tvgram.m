function [W,Wdot] = tvgram(sys,b,varargin)
% TVGRAM Controllability and observability gramians for LTV Systems on
% Finite Horizon
%
%   Wc = TVGRAM(SYS,'c') computes the controllability gramian of the
%   state-space model SYS (see TVSS) over the horizon as default time
%   scale.
%
%   Wo = TVGRAM(SYS,'o') computes its observability gramian.
%
%   In both cases, the state-space model SYS should be stable. The gramians
%   are computed by solving continuous-time matrix Lyapunov Diffrential
%   Equation:
%
%     A*Wc + Wc*A' + B*B' = -d(Wc)/dt   [Wc is controllability gramian]
%     A'*Wo + Wo*A + C'*C = -d(Wo)/dt   [Wo is observability gramian]
%
%     for continuous-time Linear-Time-Varying systems
%         dx/dt = A(t) x(t) + B(t) u(t),   y = C(t) x(t) + D(t) u(t)
% 
% For a discrete-time, time-varying model SYS:  XXX Not Implemented
%
% Example:
%
%   Time = linspace(0,5,10)';
%   AData = -5+0.1*Time.^2;
%   A = tvmat(AData,Time);
%   B = 1; C = 1; D=0;
%   G = tvss(A,B,C,D);
%   Wc = tvgram(G,'c');
%   % Options can be specified using tvodeOptions
%   tvgramopt = tvodeOptions;
%   tvgramopt.OdeOptions = odeset('RelTol',1e-5,'AbsTol',1e-8);
%   Wo = tvgram(G,'o',tvgramopt);

% Input Processing
narginchk(2,3);
if ~isa(sys,'tvss') && ~isa(sys,'tvuss')
    error('Input System must be Time Varying State Space');
end
if ~ischar(b)
    error('Second Input must be char');
end
if nargin>2
    cldeopt = varargin{1};
else
    cldeopt = [];
end

% Extract ssdata from G
[A,B,C,~] = ssdata(sys);

% Compute the gramians
switch lower(b)
    case 'c'
        % Solve clde A*Wc + Wc*A' + BB' = -d(Wc)/dt
        [W,Wdot] = clde(A,B,cldeopt);
    case 'o'
        % Solve clde A'*Wo + Wo*A + C'C = -d(Wo)/dt
        [W,Wdot] = clde(A',C',cldeopt);
    otherwise
        error('Second Input must be either ''c'' for controllability gramian or ''o'' for observability gramian.')
end
end
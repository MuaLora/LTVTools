function [Y,X] = tvlsim(G,U,x0,Opt)
% tvlsim     Simulate the response of a time-varying system
%
% [Y,X] = tvlsim(G,U) where G is a TVSS and U is a TVMAT
% [Y,X] = tvlsim(G,U,x0) where x0 is the initial condition

% XXX Simulate on time specified for U.  Check that this is 
% contained within the time horizon on which G is defined.

%% Input Processing
narginchk(2,4);
nin = nargin;
if nin==2
    Nx = order(G);
    x0 = zeros(Nx,1);
    Opt = tvodeOptions;
    Opt.OdeOptions = odeset('RelTol',1e-5,'AbsTol',1e-8);
elseif nin==3
    Opt = tvodeOptions;
    Opt.OdeOptions = odeset('RelTol',1e-5,'AbsTol',1e-8);
end
[A,B,C,D] = ssdata(G);
Tspan = [U.Time(1) U.Time(end)];

if isempty(Opt)
    Opt = tvodeOptions;
    Opt.OdeOptions = odeset('RelTol',1e-5,'AbsTol',1e-8);
elseif ~isa(Opt,'tvodeOptions') && ~(isa(Opt,'struct') &&...
        isfield(Opt,'OdeOptions') && isfield(Opt,'OdeSolver'))
    error('Options specified to tvlsim must be of type tvodeOptions.');
end

%% Simulate System (Integrate)
OdeSolver = str2func(Opt.OdeSolver);
odefh = @(t,x) LOCALderiv(t,x,A,B,U);
[t,x] = OdeSolver( odefh,Tspan,x0,Opt.OdeOptions);

% State vector X
X = tvmat(x,t);

%% Construct Outputs
U = evalt(U,t,'tvmat');
C = evalt(C,t,'tvmat');
D = evalt(D,t,'tvmat');
Y = C*X+D*U;

%% LOCALderiv
function xdot = LOCALderiv(t,x,A,B,U)
ut = evalt(U,t);    
At = evalt(A,t);
Bt = evalt(B,t);
xdot = At*x + Bt*ut;


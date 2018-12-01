function [W,Wdot] = clde(A,B,varargin)
% clde - Solve Continuous-time matrix Lyapunov Diffrential Equation
%
% [W,opt] = clde(A,B) solves the following continuous-time Lyapunov
% Matrix Diffrential Equation
%
%           A*W + W*A' + B*B' = -d(W)/dt
%           W(Tspan(1)) = F
%
% NOTE: clde uses the default horizon as specified in the A and B matrices

%% Input Processing
narginchk(2,3);
if nargin>2
    opt = varargin{1};
else
    opt = [];
end

% LIFT to tvmat
A = tvmat(A);
B = tvmat(B);

% Placeholder for checking matrix time scale and interval
% XXX - Check if A and B are with the same time horizon
% XXX - Check for the need of interpolation
if ~isequal(A.Time(end),B.Time(end)) || ~isequal(diff(A.Time(end)),diff(B.Time(end)))
    error('The input matrices must be of the same time horizon and time intervals');
end

% Specify Tspan (Backwards in Time)
Tf = A.Time(end);
Tspan = [Tf 0];

%% Set Default Values
% Initialize W(Tf) as F being all zeros
Nx = size(A,1);
F = zeros(Nx);
F = F(:);

% PARSE opt to OdeSolver and OdeOptions
if isa(opt,'tvodeOptions') || (isa(opt,'struct') &&...
        isfield(opt,'OdeOptions') && isfield(opt,'OdeSolver'))
    OdeOpt = opt.OdeOptions;
    OdeSolverStr = opt.OdeSolver;
elseif isempty(opt) % If isempty(opt) then specify default values for OdeSolver and
    opt = tvodeOptions;
    OdeOpt = odeset('RelTol',1e-5,'AbsTol',1e-8);
    OdeSolverStr = opt.OdeSolver; % Default
else
    error('The options must be of type tvodeOptions or a structure with OdeOptions and OdeSolver as fields.');
end

% XXX - The code below overwrites any Event specified by the user.
OdeOpt = odeset(OdeOpt,'Events',@LOCALevents);
OdeSolver = str2func(OdeSolverStr);

%% Create function handles for ODE solver
Afh = LOCALtv2fh(A);
Bfh = LOCALtv2fh(B);
odefh = @(t,W) LOCALWdot(t,W,Afh,Bfh);

%% Solve Lyapunov Diffrential Equation
% Note: Setting the time span as [T, 0] indicates that OdeSolver should
% integrate backwards from the boundary condition P(T)=F.
% warning off;
%[t,P] = ode45( odefh,tspan,F,odeopts);
warning('off',['MATLAB:' OdeSolverStr ':IntegrationTolNotMet']);
sol = OdeSolver(odefh,Tspan,F,OdeOpt);
warning('on',['MATLAB:' OdeSolverStr ':IntegrationTolNotMet']);

%% Construct Outputs
% Note: sol = ode45() returns a solution structure that can be
% evaluated with deval. The alternative syntax [t,P]=ode45() refines
% the solution time grid by a factor of 4.  In other words,
%   t(1)=sol.x(1), t(5)=sol.x(2), t(9)=sol.x(3).
% The documentation states this default refinement is used because
% ODE45 can take large time steps.

% Time Vector
if numel(Tspan)==2
    if isempty(OdeOpt.Refine)
        Nr = 4;    % This is the default ODESET option
    else
        Nr = OdeOpt.Refine;
    end
    fac = ((1:Nr)-1)/Nr;
    tSol = sort(sol.x);
    Nsol = numel(tSol);
    tDiff = diff(tSol);
    t = zeros(Nr,Nsol-1);
    t(1,:) = tSol(1:end-1);
    for i=2:Nr
        t(i,:) = tSol(1:end-1) + fac(i)*tDiff;
    end
    t = [t(:); tSol(end)];
    Nt = numel(t);
else
    t = sort(Tspan);
end

% Compute W and Wdot
nout = nargout;
if isequal(nout,1)
    W = deval(sol,t);
    W = tvmat( reshape(W,[Nx Nx Nt]) , t );
else
    % Two options for computing Pdot:
    %  A) DEVAL computes Pdot as "the first derivative of the polynomial
    %        interpolating the solution"
    %  B) Directly call the ODEFH to compute Pdot given P
    % Option B seems like it should be more accurate (but perhaps
    % more computationally costly?)
    
    [W,Wdot] = deval(sol,t);
    W = tvmat( reshape(W,[Nx Nx Nt]) , t );
    Wdot = tvmat( reshape(Wdot,[Nx Nx Nt]) , t );
end
end

%% LOCAL Function
function Mfh = LOCALtv2fh(M)

if isempty(M)
    % Mfh = M;  % Defines Mfh as an empty TVMAT which is slower than
    Mfh = [];
elseif M.isConstant
    Mfh = @(t) M.Data;
else
    if isequal(M.InterpolationMethod,'Spline') && isempty(M.SplineData)
        M = getSplineData(M);
    end
    Mfh = @(t) evalt(M,t);
end
end

%% LOCAL Function
% This can be improved by:
%  1) Accounting for the symmetry in P;
%  2) Possibly doing some initial transformation on (A,B,Q,S,R)
%     to put the data in a form for more efficient integration.
%     [Perhaps look at work by Laub, others for AREs]
%  3) Connection to integrating Hamiltonian?  Look at the PLTV literature.
function Wdot = LOCALWdot(t,W,Afh,Bfh)

% Get Lyapunov state matrices
A = Afh(t);
B = Bfh(t);

% Convert W from column to matrix
Nx = size(A,1);
W = reshape(W,[Nx, Nx]);

% Compute Wdot (Lyapunov Diffrential Equation)
Wdot = - A*W - W*A' - B*B';

% Convert Pdot from matrix back to column
Wdot = Wdot(:);
end

%% LOCAL Function
function [value,isterminal,direction] = LOCALevents(~,W)

if any(isnan(W(:)))
    value = 0;
else
    value = 1;
end
isterminal = 1;
direction = 0;
end

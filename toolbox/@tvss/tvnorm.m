function [g,d,info] = tvnorm(G,Tspan,varargin)
%% tvnorm Computes the time varying system norm
%
% [g,d,info] = tvnorm(G,Tspan,NormType,Opt)
%
% G is a time varying state space object
% NormType = 'L2toL2' or 'L2toE'. This input is optional with default 'L2toL2'.

%% Input Processing

% Check # of inputs
narginchk(2,4);
nin = nargin;
nout = nargout;

% Parse inputs
NormType = [];
Opt = [];
if nin==3
    if isa(varargin{end},'tvnormOptions')
        Opt = varargin{1};
    else
        NormType = varargin{1};
    end
elseif nin==4
    [NormType,Opt] = deal(varargin{:});
end

% Default Inputs
if isempty(NormType)
    NormType = 'L2toL2';
end
if isempty(Opt)
    Opt = tvnormOptions;
end

% Process Time Span
if isscalar(Tspan)
    T0=0;
    Tf=Tspan;
else
    T0 = Tspan(1);
    Tf = Tspan(2);
end

% Process Options
if isequal(Opt.Display,'on')
    DispFlag = true;
else
    DispFlag = false;
end

% Update crdeOpt based on tvnormoptions default
crdeOpt.OdeOptions = Opt.OdeOptions;
crdeOpt.OdeSolver = Opt.OdeSolver;

%% Build Cost Function Matrices
% Note:  R(t,g) = R0(t) - g^2 R1(t)
[A,B,C,D] = ssdata(G);
[Nx,Nd] = size(B);
E = [];
R1 = eye(Nd);

switch NormType
    case 'L2toL2'
        Q = C'*C;
        S = C'*D;
        R0 = D'*D;
        F = zeros(Nx);
        
    case 'L2toE'
        % XXX Verify D(t)==0 for L2 to E
        % XXX This set-up requires R(t) to be inverted at each time step.
        % However, no inversion is actually required since we can set
        % R=I and B'=B/gamma.  This should speed up the integration.
        Q = zeros(Nx);
        S = [];
        R0 = zeros(Nd);
        CTf = evalt(C,Tf);
        F = CTf'*CTf;
        
    otherwise
        error('NormType must be set to either L2toL2 or L2toE.');
end

%% Lower Bound Phase
%  Require R(t,g)=R0(t)-g^2*R1(t)<0 for all t in [0,T]
% XXX This assumes the time grid is sufficiently fine.
if isequal(NormType,'L2toL2')
    ev = eig(R0);
    evmax = max( ev.Data(:) );  % XXX Implement a TVMAX function?
    gLow = sqrt(evmax);
else
    gLow = 0;
end
gLow = max(gLow,Opt.Bounds(1));
if DispFlag
    fprintf('\n Lower Bound = %4.3f',gLow);
end

%% Upper Bound Phase
PUpp = []; PdotUpp = []; solUpp = [];
PLow = []; PdotLow = []; solLow = [];
if isfinite(Opt.Bounds(2))
    % User specified a (finite) upper bound.
    % Verify (or disprove) this upper bound.
    gTry = Opt.Bounds(2);
    R = R0-gTry^2*R1;
    if nout==1
        P = crde(A,B,Q,R,S,E,F,[Tf T0],crdeOpt);
        Pdot = []; sol = [];
    else
        [P,~,Pdot,sol] = crde(A,B,Q,R,S,E,F,[Tf T0],crdeOpt);
    end
    
    % Check convergence of P
    if P.Time(1)>T0
        % P did not converge
        gUpp = inf;
        gLow = gTry;
        PLow = P; PdotLow = Pdot; solLow = sol;
    else
        % P converged
        gUpp = gTry;
        PUpp = P; PdotUpp = Pdot; solUpp = sol;
    end
else
    % User did not specify a (finite) upper bound.
    % Attempt to determine a finite upper bound on performance.
    gFac = 10;
    gUpp = gFac*gLow+1;  % XXX Better choice?
    
    cnt = 0;
    cntMax = 8;
    haveUpper = false;
    while ~haveUpper && cnt<cntMax
        % Pick gamma
        cnt = cnt+1;
        gTry = gUpp;
        
        % Solve LTV Riccati Equation
        R = R0-gTry^2*R1;
        if nout==1
            P = crde(A,B,Q,R,S,E,F,[Tf T0],crdeOpt);
            Pdot = []; sol = [];
        else
            [P,~,Pdot,sol] = crde(A,B,Q,R,S,E,F,[Tf T0],crdeOpt);
        end
        
        % Check convergence of P
        if P.Time(1)>T0
            % P did not converge
            gUpp = gFac*gUpp;
            gLow = gTry;
            PLow = P; PdotLow = Pdot; solLow = sol;
        else
            % P converged
            haveUpper = true;
            gUpp = gTry;
            PUpp = P; PdotUpp = Pdot; solUpp = sol;
        end
    end
    
    if ~haveUpper
        % Could not find a finite upper bound
        if DispFlag
            fprintf(['\n Could not find a finite upper bound.' ...
                ' Infeasible at gTry = %4.1f\n'],gUpp);
        end
        gUpp = inf;
    else
        if DispFlag
            fprintf('\n Lower Bound = %4.3f \t Upper Bound = %4.3f',...
                gLow,gUpp);
        end
    end
end

%% Bisection Phase
AbsTol = Opt.AbsTol;
RelTol = Opt.RelTol;
if isfinite(gUpp)
    if DispFlag
        fprintf('\n');
    end
    while (gUpp - gLow > RelTol*gUpp+AbsTol)
        gTry = (gUpp+gLow)/2;
        
        % Solve LTV Riccati Equation
        R = R0-gTry^2*eye(Nd);
        if nout==1
            P = crde(A,B,Q,R,S,E,F,[Tf T0],crdeOpt);
            Pdot = []; sol = [];
        else
            [P,~,Pdot,sol] = crde(A,B,Q,R,S,E,F,[Tf T0],crdeOpt);
        end
        
        % Check convergence of P
        if P.Time(1)>T0
            % P did not converge
            gLow = gTry;
            PLow = P; PdotLow = Pdot; solLow = sol;
            if DispFlag
                fprintf(' gTry = %4.3f \t tConv = %4.3f (N) \t gLow = %4.3f \t gUpp = %4.3f\n', ...
                    gTry,P.Time(1),gLow,gUpp);
            end
        else
            % P converged
            gUpp = gTry;
            PUpp = P; PdotUpp = Pdot; solUpp = sol;
            if DispFlag
                fprintf(' gTry = %4.3f \t tConv= %4.3f  (Y) \t gLow = %4.3f \t gUpp = %4.3f\n', ...
                    gTry,P.Time(1),gLow,gUpp);
            end
        end
        
    end
end

%% Store Final Result
g = [gUpp, gLow];

info.Lower.Gain = gLow;
info.Lower.P = PLow;
info.Lower.Pdot = PdotLow;
info.Lower.sol = solLow;

info.Upper.Gain = gUpp;
info.Upper.P = PUpp;
info.Upper.Pdot = PdotUpp;
info.Upper.sol = solUpp;


%% Construct Worst-Case Input
% The construction is based on the two-point boundary value problem
% (TPBVP) related to the (Q,S,R,F) cost.  This particular implementation
% using the transformed Hamiltonian dynamics was developed by A. Iannelli
% and P. Seiler in September 2017.
d = [];
if nout>=2 && ~isempty(PLow)
    % Compute largest eigenvalue of PLow
    T = PLow.Time;
    [evecP,evalP] = eig( evalt(PLow,T(1)) );
    [emax,idx] = max( diag(evalP) );
    vmax = evecP(:,idx);
    
    %norm( evalt(PLow,T(1))*vmax - emax*vmax )  %norm(emax*vmax)]
    
    
    % Evaluate State/Cost Matrices on Same Time Grid as P
    R = tvmat( R0-gLow^2*R1 );
    A = evalt(A,T,'tvmat');
    B = evalt(B,T,'tvmat');
    if isequal(NormType,'L2toL2')
        R = evalt(R,T,'tvmat');
        S = evalt(tvmat(S),T,'tvmat');
    else
        S = tvmat(zeros(Nx,Nd));
    end
    
    % Simulate Transformed Hamiltonian system
    M = inv(R)*(PLow*B+S)';
    H11 = A-B*M;
    odefh = @(t,E) evalt(H11,t)*E;
    E0 = vmax/emax;
    [tE,E] = ode45( odefh,T,E0);
    E = tvmat(E,tE);
    
    % Construct disturbance
    % Note: This calculation uses a non-convergent Riccati solution and
    % hence the disturbance starts at PLow.Time(1) > T0.
    d = -M*E;
    
    % Append zero input to disturbance on the window [PLow.Time(1) T0)
    dTime = d.Time;
    dData = d.Data;
    
    dTime = [T0; T0+0.999*(dTime(1)-T0); dTime];
    dData = cat(3,zeros(Nd,1,2), dData);
    d = tvmat(dData,dTime);
end
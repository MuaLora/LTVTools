function H = evalt(G,T,tvflag)
% EVALT Evaluate a TVSS at specified time(s)
%
% H = EVALT(G,T) evaluates a TVSS G on the vector of times T using the
% interpolation method specified by G.InterpolationMethod.  H is returned
% as an SS array.
%
% H = EVALT(G,T,'TVSS') evaluates a TVSS G on the vector of times T.
% H is returned as a TVSS object.

%% Default is to return a double unless user specifies "TVSS" option
if nargin==3 && strcmpi(tvflag,'TVSS')
    tvflag = true;
    tvstr = {'TVMAT'};
else
    tvflag = false;
    tvstr = {};
end

%% Evaluate State Matrices
[A,B,C,D] = ssdata(G);
At = evalt(A,T,tvstr{:});
Bt = evalt(B,T,tvstr{:});
Ct = evalt(C,T,tvstr{:});
Dt = evalt(D,T,tvstr{:});

%% Package up result
if tvflag
    H = tvss(At,Bt,Ct,Dt);
else
    H = ss(At,Bt,Ct,Dt);
end
function [gFinal,wcinfo] = tvrobL2toL2(G,v,p,Tf,tlmi,tSp)
%% tvrobL2toL2  Worst-case L2-to-L2 norm
% [g,wcinfo] = tvrobL2toE(G,v,p,Tf,tlmi,tSp)
%
% Inputs
%   G - Nominal LTV system (as a TVSS). The uncertain system is
%       Fu(Delta,G) where Delta is a SISO, LTI, uncertainty
%       satisyfing the unit norm bound ||Delta|| <= 1.
%   (v,p) - The IQC multiplier is blkdiag(Psiv'*X11*Psiv,-Psiv'*X11*Psiv)
%           where Psiv:=[1; 1/(s-p); ...; 1/(s-p)^v] and X11>=0.
%   Tf - Analysis is performed on the horizon [0,Tf]
%   tlmi - Nlmi-by-1 vector of time points to enforce the LMI
%   tSp - Nsp-by-1 vector of time points for the spline basis
%
% Outputs
%   g - Gain Upper bound
%   wcinfo - Structure of solution info

%% Input processing
t0 = tic;
tlmi = tlmi(:);
DispFlag = true;

%% Create Spline Basis Functions
Nsp = numel(tSp);
Ps = [];
for i=1:Nsp
    ei = zeros(Nsp,1);
    ei(i) = 1;
    Ps = [Ps; tvmat(ei,tSp,'Spline')];
end

%% Iteration
Niter = 10;         % Maximum # of iterations
glmi = zeros(Niter,1);
grde = zeros(Niter,1);
for i=1:Niter
    if DispFlag
        fprintf('\n Iteration # = %d',i);
    end
    
    % Evaluate data on LMI time grid
    G1 = evalt(G,tlmi,'tvss');
    Ps1 = evalt(Ps,tlmi,'tvmat');
    Psdot1 = tvdiff(Ps,tlmi);
    if i>1 && isfinite(grde(i-1))
        % Use RDE solution from previous iteration (if grde<inf)
        [Pm,Pmdot] = evalRDE(solrdePrev,tlmi);
    else
        Pm = tvmat; Pmdot = Pm;
    end
        
    % Finite Horizon: LMI Condition
    [glmi(i),X11,Plmi,Pdotlmi,Y,lamP] = ...
        FHL2toL2lmi(G1,v,p,Ps1,Psdot1,Pm,Pmdot);
    if DispFlag
        fprintf('\n LMI Gain Bound = %4.4f,',glmi(i));
    end
    LMIinfo = struct('X11',X11,'P',Plmi,'Pdot',Pdotlmi,'Y',Y,'lamP',lamP);
    
    % Evaluate LMI solution (for debugging) to show eLMI1<=0
    % [LMI1,eLMI1] = evalL2toL2LMI(G1,v,p,glmi(i),X11,Plmi,Pdotlmi);
    % figure; plot(eLMI1.Time,eLMI1.Data(:));
    
    % Finite Horizon: RDE + Bisection
    [gbnds,RDEinfo] = FHL2toL2rde(G,v,p,Tf,X11);
    grde(i) = gbnds(1);
    if isfinite(grde(i))
        trde = fliplr(RDEinfo.Upper.sol.x);
        solrde = RDEinfo.Upper.sol;
    else
        trde = fliplr(RDEinfo.Lower.sol.x);
        solrde = RDEinfo.Lower.sol;
    end
    
    if DispFlag
        fprintf('\t RDE Gain Bound = %4.4f,',grde(i));
        if isinf(grde(i))
            fprintf('\t trde = %4.3f\n',trde(1));
        else
            fprintf('\n');
        end
    end    

    % Termination Condition: LMI and RDE gains are "close"
    if abs(grde(i) - glmi(i)) < 5e-3*glmi(i)
        break;
    end
    
    if glmi(i) < grde(i)
        % LMI cost was better indicating the LMI time grid was too coarse
        
        % Create dense time grid
        tDense = union(trde,tlmi);
        t1 = trde(1);
        if t1>0
            % RDE solution diverged at t1. Add time points in [0,t1]
            tdensity = numel( tDense(tDense>=t1) ) / (Tf-t1);
            tadd = linspace(0,t1, ceil(t1*tdensity) );
            tDense = union(tDense,tadd);
        end
        
        % Evaluate LMI using (P,Pdot,X11,g) but on dense time grid
        G2 = evalt(G,tDense,'tvss');
        Ps2 = evalt(Ps,tDense,'tvmat');
        Psdot2 = tvdiff(Ps,tDense);        
        if i>1 && isfinite(grde(i-1))
            % Convergent RDE solution on previous iteration
            [Pm2,Pmdot2] = evalRDE(solrdePrev,tDense);
        else
            % No convergent RDE solution on previous iteration
            Pm2 = tvmat; Pmdot2 = Pm2;
        end
        
        [Pval,Pdotval] = evalP(Ps2,Psdot2,Pm2,Pmdot2,Y,lamP);
        [~,eLMI2] = evalL2toL2lmi(G2,v,p,glmi(i),X11,Pval,Pdotval);
        
        % Update the LMI time grid.
        eLMI2 = eLMI2.Data(:);
        emax = max(eLMI2);
        addidx = [];
        for j=1:numel(tlmi)-1
            idx = find( tDense>tlmi(j) & tDense<tlmi(j+1) );
            if ~isempty(idx)
                [~,idx2] = max( eLMI2(idx) );
                if eLMI2(idx(idx2))> 0.25*emax
                    addidx = [addidx; idx(idx2)];
                end
            end
        end
        tadd = tDense(addidx);
        
        if DispFlag
            fprintf(' Adding t = ');
            for j=1:numel(tadd)
                fprintf('%4.3f, ',tadd(j));
            end
            fprintf('\n');
        end
        
%         figure(3)
%         plot(tlmi,eLMI1,'b',tlmi,eLMI1,'bx',tt,eLMI2,'r',...
%             tadd,eLMI2(addidx),'rx')
%         ylim(max(eLMI2(:))*[-1 1]);
%         drawnow;
        
        tlmiNew = sort([tlmi; tadd]);
        if isequal(tlmiNew,tlmi)
            % Terminating - Same Time Grid
            break;
        end
        tlmi = tlmiNew;
    end
    
    % Store Iteration Info
    solrdePrev = solrde;
    wcinfo(i) = struct('TotalTime',toc(t0),'LMI',LMIinfo,'RDE',RDEinfo);
end

% Store Iteration Info (if iteration terminated before Niter)
if i<Niter
    wcinfo(i) = struct('TotalTime',toc(t0),'LMI',LMIinfo,'RDE',RDEinfo);
end

% Store outputs
gFinal = min(grde(1:i)); % grde(i);
tL2toE = toc(t0);
if DispFlag
    fprintf('\n Final Results:');
    fprintf(' Robust gain = %4.4f',gFinal);
    fprintf('  Total Comp Time = %4.4f\n',tL2toE);
end


%% LOCAL Function: evalP
function [P,Pdot] = evalP(Ps,Psdot,Pm,Pmdot,Y,lamP)
% Combine scalar and matrix bases functions
%    P = ps1*Y1+ ... psN*YN + lamP*Pm
% The derivative Pdot is computed simiarly.

% Sum terms with scalar basis functions
Ps = Ps(:);
Ns = size(Ps,1);
P = Ps(1)*Y(:,:,1);
Pdot = Psdot(1)*Y(:,:,1);
for j=2:Ns
    P = P + Ps(j)*Y(:,:,j);
    Pdot = Pdot + Psdot(j)*Y(:,:,j);
end

% Include matrix basis function (if it exists)
if ~isempty(Pm)
    % XXX
    Pm.InterpolationMethod = 'Spline';
    Pmdot.InterpolationMethod = 'Spline';
    
    P = P + lamP*Pm;
    Pdot = Pdot + lamP*Pmdot;
end


%% LOCAL Function: evalRDE
function [P,Pdot] = evalRDE(solrde,t)
% Evaluate RDE solution on a specified time grid.

% Dimensions
Nt = numel(t);
Nx = round(sqrt(  size(solrde.y,1) ));

% Using deval to evaluate Pdot is different (less accurate?) than
% directly calling the ODEFH to compute the derivative (as done below)
%     [P,Pdot] = deval(solrde,t);
%     P = reshape( P, [Nx Nx Nt]);
%     Pdot = reshape( Pdot, [Nx Nx Nt]);

% Construct P
P = deval(solrde,t);
P = reshape(P,[Nx Nx Nt]);

% Evaluate Pdot
odefh = solrde.extdata.odefun;
Pdot = zeros(Nx,Nx,Nt);
for i=1:Nt
    Pdot(:,:,i)=reshape( odefh(t(i),P(:,:,i)), [Nx Nx]);
end

% Ensure symmetry
for i=1:Nt
    P(:,:,i) = ( P(:,:,i)+P(:,:,i)' )/2;
    Pdot(:,:,i) = ( Pdot(:,:,i)+Pdot(:,:,i)' )/2;
end

% Return as TVMAT
P = tvmat(P,t);
Pdot = tvmat(Pdot,t);


%% LOCAL Function: evalL2toL2lmi
function [LMI,eLMI] = evalL2toL2lmi(G,v,p,g,X11,P,Pdot)
% This function evaluates the finite-horizon induced L2 gain condition
% on a time grid for given values of the IQC, storage function, and gain.
% **The function assumes G is LTV and Delta is a SISO, unit norm bounded,
% LTI uncertainty.
%
% NOTE - All data (G,P,Pdot) must have the same time grid t

% Form Extended System
[AAe,BBe,CCe,DDe,Nz,~,Nd] = ExtSystem(G,v,p);
n = size(AAe,1);
CC1e = CCe(1:Nz,:); 
DD1e = DDe(1:Nz,:);
CC2e = CCe(Nz+1:end,:); 
DD2e = DDe(Nz+1:end,:);  

% Evaluate LMI
P.InterpolationMethod = 'Linear';
Pdot.InterpolationMethod = 'Linear';

X = blkdiag(X11,-X11);
gsq = g^2;
LMI = [AAe'*P+P*AAe+Pdot P*BBe; BBe'*P zeros(Nd+1)] ...
    + blkdiag( zeros(n+1), -gsq*eye(Nd) ) + [CC1e DD1e]'*X*[CC1e DD1e] ...
    + [CC2e DD2e]'*[CC2e DD2e];

% Evaluate max real eigenvalue
eLMI = max(real( eig(LMI) ));


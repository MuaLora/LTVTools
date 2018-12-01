function B = evalt(A,T,tvflag)
% EVALT Evaluate a TVMAT at specified time(s)
%
% B = EVALT(A,T) evaluates a TVMAT A on the vector of times T using the
% interpolation method specified by A.InterpolationMethod.  B is returned
% as a double array.
%
% B = EVALT(A,T,'TVMAT') evaluates a TVMAT A on the vector of times T.
% B is returned as a TVMAT object.

%% Default is to return a double unless user specifies "TVMAT" option
if nargin==3 && strcmpi(tvflag,'TVMAT')
    tvflag = true;
else
    tvflag = false;
end

%% Handle Constant Case
AData = A.Data;
Nd = ndims(AData);
Nt = numel(T);
if A.isConstant
    BData = repmat(AData,[ones(1,Nd) Nt]);
    if tvflag
        B = tvmat(BData,T);
    else
        B = BData;
    end
    return
end

%% Handle Varying Case
ATime = A.Time;
AIM = A.InterpolationMethod;

% Pre-process TVMAT Data
ADiff = diff(ATime);
NATime = numel(ATime);
nad = ndims(AData)-3;
if nad==0
    id =  cell(1,0);
else
    id = repmat({':'},1,nad);
end

% Perform Interpolation
switch AIM
    case 'Linear'
        for i=1:Nt
            [k,alpha] = LOCALfindslotalpha(NATime,ATime,T(i),ADiff);
            if alpha==0
                Bi = AData(:,:,id{:},k);
            else
                Bi = (1-alpha)*AData(:,:,id{:},k) + ...
                    alpha*AData(:,:,id{:},k+1);
            end
            if i==1
                BData = zeros( [numel(Bi) Nt] );
            end
            BData(:,i) = Bi(:);
        end
        BData = reshape(BData,[size(Bi) Nt]);
    case 'Flat'
        for i=1:Nt
            k = LOCALfindslotalpha(NATime,ATime,T(i),ADiff);
            Bi = AData(:,:,id{:},k);
            if i==1
                BData = zeros( [numel(Bi) Nt] );
            end
            BData(:,i) = Bi(:);
        end
        BData = reshape(BData,[size(Bi) Nt]);
        
    case 'Nearest'
        for i=1:Nt
            [k,alpha] = LOCALfindslotalpha(NATime,ATime,T(i),ADiff);
            if alpha>=0.5
                k=k+1;
            end
            Bi = AData(:,:,id{:},k);
            if i==1
                BData = zeros( [numel(Bi) Nt] );
            end
            BData(:,i) = Bi(:);
        end
        BData = reshape(BData,[size(Bi) Nt]);
        
    case 'Spline'
        % Currently SplineData is empty and it computed in every call
        % to EVALT. A syntax to precompute Spline Data is:
        %    A = getSplineData(A);
        % This allows for more efficient, repeated calls to EVALT.
        ASplineData = A.SplineData;
        if isempty(ASplineData)
            ASplineData = getSplineData(AData,ATime);
        end
        
        % XXX The code below allows time data outside Time range but
        % clips to nearest time point. Error out instead?
        for i=1:Nt
            if T(i)<=ATime(1)
                k=1;
                dt = 0;
            elseif T(i)>=ATime(end)
                k=numel(ATime)-1;
                dt = ATime(end)-ATime(end-1);
            else
                k = find( T(i)>=ATime(1:end-1) & T(i)<ATime(2:end) );
                k = k(1);
                dt = T(i)-ATime(k);
            end
            aidx = 4*(k-1)+(1:4);
            if dt==0
                Bi = ASplineData(aidx(1),:);
            else
                Bi = [1 dt dt^2 dt^3]*ASplineData(aidx,:);
            end
            
            if i==1
                BData = zeros( [numel(Bi) Nt] );
            end
            BData(:,i) = Bi(:);
        end
        BData = reshape(BData,[size(AData(:,:,id{:},1)) Nt]);
end

%% Convert output to a TVMAT if option specified
if tvflag
    B = tvmat(BData,T,AIM);
else
    B = BData;
end



function [k,alpha] = LOCALfindslotalpha(N,vec,val,dvec)
% N integer
% vec 1-by-N (or N-by-1), sorted
% val 1-by-1
% dvec = diff(vec)

% XXX This allows time data outside range but clips to nearest time point.
% Error out instead?

k = max(find(val>=vec));   %#ok<MXFND> % don't follow advice - it is slower.
if ~isempty(k)
    if k<N
        alpha = (val - vec(k))/dvec(k);
    else
        alpha = 0;
    end
else
    k = 1;
    alpha = 0;
end



% TVMAT Create a time-varying matrix
%
%   M = tvmat(Data,Time) creates a time-varying matrix. Time is an Nt-by-1
%   vector of increasing time values.  If Data is 1-by-Nt then M is a
%   time-varying scalar with value Data(i) at Time(i). If Data is N1-by-Nt
%   then M is an N1-by-1 time-varying vector with value Data(:,i) at
%   Time(i). In general, if Data is an N1-by-N2-by-...-by-Nad-by-Nt array
%   then M is an N1-by-N2-by-...-by-Nad time-varying array with value
%   Data(:,...,:,i) at Time(i).
%
%   M = tvmat(Data,Time,InterpolationMethod) specifies the method to
%   interpolate the data between the specified time grid points.  The
%   method can be 'Linear', 'Flat', 'Nearest', or 'Spline'.
%   Default is InterpolationMethod = 'Linear'.
%
%   M = tvmat(Data) generates a constant matrix.
%
%   % EXAMPLE: (CUT/PASTE)
%   % Create a 2-by-2 matrix defined on a 1-dimensional grid
%   M1 = tvmat(5);
%   Time = linspace(0,5,10)';
%   Data = randn(1,1,10);
%   M2 = tvmat(Data,Time);
%   M3 = M1+M2;
%   plot(M2.Time,M2,'b',M3.Time,M3,'r--');
%
% See also:

classdef (InferiorClasses={?frd, ?ss,?tf,?zpk,?ureal,?ucomplex,?ucomplexm,...
        ?ultidyn,?udyn,?umat,?uss,?ufrd}) tvmat
    % Class definition for time-varying matrix
    
    properties
        Data = [];
        Time = [-inf,inf];
        InterpolationMethod = 'Linear';  % Flat, Nearest, Spline,
        % The method names align (except for spline) with those used in
        % the LPV Simulink block.
        %
        % XXX
        % Add flag to compute spline data?  (We don't want to compute
        %   this every time the constructor is called as this causes
        %   all operations to be slow.  However, we also need some
        %   method to say that we want to compute it otherwise it will
        %   be computed everytime EVALT is called)
        % Add TimeUnits,
        % Store offset data? (e.g. as returned by linearize)
    end
    
    properties (SetAccess = immutable)
        isConstant = false;
    end
    
    properties (SetAccess = private)
        % Note: This class has a subsasgn method which, in general, by-
        % passes "SetAccess = private". The private SetAccess is actually
        % enforced by logic in subsasgn.
        SplineData = [];
    end
    
    methods
        % Constructor
        % XXX Add undocumented flag to skip error checking?
        function obj = tvmat(Data,Time,InterpolationMethod)
            nin = nargin;
            if nin==0
                % Return default object
                return
            elseif nin==1
                if isa(Data,'tvmat')
                    obj = Data;
                    return
                else
                    obj.Data = Data;
                    obj.isConstant = true;
                end
            elseif nin==2
                obj.Data = Data;
                obj.Time = Time(:);
            elseif nin==3
                obj.Data = Data;
                obj.Time = Time(:);
                obj.InterpolationMethod = InterpolationMethod;
            else
                error('Invalid syntax for the "tvmat" command.');
            end
            
            Nt = numel(obj.Time);
            if ~obj.isConstant
                % Code below allows Data to be Nt-by-Nd or Nd-by-Nt.
                % It reshapes Data to Nd-by-1-by-Nt.
                if ismatrix(Data)
                    szData = size(Data);
                    if szData(2)==Nt
                        obj.Data = reshape(Data,[szData(1) 1 Nt]);
                    elseif szData(1)==Nt
                        obj.Data = reshape(Data',[szData(2) 1 Nt]);
                    end
                end
            end
            
            % Squeeze extra singleton dimensions
            if ~obj.isConstant && ~all(size(Data)==0)
                szData = size(obj.Data);
                szData(end) = 1;
                idx = find( diff(szData)~= 0 );
                Nd = max(2, numel(idx) );
                obj.Data = reshape(obj.Data, [szData(1:Nd) Nt]);
            end
            
            % Use isvalid to perform all error checking
            [pflag,errstr] = isvalid(obj);
            if pflag==0
                error(errstr);
            end
        end
        
        % isvalid: Determine if TVMAT object has valid data
        function [pflag,errstr] = isvalid(obj)
            % ISVALID Determine if TVMAT object is valid.
            errstr = [];
            pflag = 1;
            
            % Check Data Class
            Data = obj.Data;
            if ~( isa(Data,'double') || isa(Data,'logical')  )
                pflag = 0;
                errstr = 'Data must be a double or logical.';
                return
            end
            
            % Check Time
            Time = obj.Time;
            Nt = numel(Time);
            if ~isa(Time,'double') || Nt<2 || any( diff(Time)< 0 )
                pflag = 0;
                errstr = ['Time vector must be non-decreasing vector ' ...
                    'of doubles with length>=2.'];
                return
            end
            
            % Check Time/Data Consistency
            szData = size(obj.Data);
            if obj.isConstant && ~all(Time==[-inf inf])
                pflag = 0;
                errstr = ['Time vector must be [-inf inf] for constant TVMATs.'];
            end
            if ~obj.isConstant && ~all(size(Data)==0) && ...
                    ~(numel(szData)>2 && szData(end)==Nt )
                pflag = 0;
                errstr = 'Dimensions of Time and Data are incompatible.';
            end
            
            % Check IM
            IM = obj.InterpolationMethod;
            if ~any( strcmpi(IM,{'Flat'; 'Nearest';'Linear'; 'Spline'}) )
                pflag = 0;
                errstr = ['InterpolationMethod must be "Flat", '...
                    '"Nearest", "Linear", or "Spline".'];
            end
            if pflag==0 && nargout==0
                error(errstr);
            end
        end
        
        % % XXX Create Display
        % function s=display(obj)
        % end
        %
        % XXX numel, numArgumentsFromSubscript, get , set
        
        function out = double(obj)
            out = obj.Data;
        end
        
        function out = end(obj,slot,nslots)
            sz = size(obj);
            if slot < nslots
                out = sz(slot);
            else
                out = prod(sz(slot:end));
            end
        end
        
        function out = length(obj)
            out = size(obj);
            out = max(out);
        end
        
        function out = isempty(obj)
            %szo = size(obj); out = any(~szo);
            out = isempty(obj.Data);
        end
        
        function out = iscolumn(obj)
            out = size(obj,2)==1;
        end
        
        function out = isrow(obj)
            out = size(obj,1)==1;
        end
        
        function out = isscalar(obj)
            szo = size(obj);
            out = (szo(1)==1 && szo(2)==1);
        end
        
        function out = isvector(obj)
            szo = size(obj);
            out = (szo(1)==1 || szo(2)==1);
        end
        
        function out = ndims(obj)
            if obj.isConstant
                out = ndims( obj.Data );
            else
                out = ndims( obj.Data ) - 1;
            end
        end
        
        function varargout = size(obj,varargin)
            nout = nargout;
            if nout==0
                varargout = cell(1,1);
            else
                varargout = cell(1,nout);
            end
            if obj.isConstant
                [varargout{:}] = size( obj.Data, varargin{:});
            else
                Data = obj.Data;
                nd = ndims(Data)-3;
                id = repmat({':'},1,nd);
                [varargout{:}] = size( Data(:,:,id{:},1), varargin{:});
            end
        end
        
        
        %--------------Unary Operations: Element-by-Element
        function out = abs(obj)
            out = tvunopebe(@abs,obj);
        end
        function out = acos(obj)
            out = tvunopebe(@acos,obj);
        end
        function out = acosd(obj)
            out = tvunopebe(@acosd,obj);
        end
        function out = acosh(obj)
            out = tvunopebe(@acosh,obj);
        end
        function out = acot(obj)
            out = tvunopebe(@acot,obj);
        end
        function out = acotd(obj)
            out = tvunopebe(@acotd,obj);
        end
        function out = acoth(obj)
            out = tvunopebe(@acoth,obj);
        end
        function out = acsc(obj)
            out = tvunopebe(@acsc,obj);
        end
        function out = acscd(obj)
            out = tvunopebe(@acscd,obj);
        end
        function out = acsch(obj)
            out = tvunopebe(@acsch,obj);
        end
        function out = angle(obj)
            out = tvunopebe(@angle,obj);
        end
        function out = asec(obj)
            out = tvunopebe(@asec,obj);
        end
        function out = asecd(obj)
            out = tvunopebe(@asecd,obj);
        end
        function out = asech(obj)
            out = tvunopebe(@asech,obj);
        end
        function out = asin(obj)
            out = tvunopebe(@asin,obj);
        end
        function out = asind(obj)
            out = tvunopebe(@asind,obj);
        end
        function out = asinh(obj)
            out = tvunopebe(@asinh,obj);
        end
        function out = atan(obj)
            out = tvunopebe(@atan,obj);
        end
        function out = atand(obj)
            out = tvunopebe(@atand,obj);
        end
        function out = atanh(obj)
            out = tvunopebe(@atanh,obj);
        end
        function out = ceil(obj)
            out = tvunopebe(@ceil,obj);
        end
        function out = conj(obj)
            out = tvunopebe(@conj,obj);
        end
        function out = cos(obj)
            out = tvunopebe(@cos,obj);
        end
        function out = cosd(obj)
            out = tvunopebe(@cosd,obj);
        end
        function out = cosh(obj)
            out = tvunopebe(@cosh,obj);
        end
        function out = cot(obj)
            out = tvunopebe(@cot,obj);
        end
        function out = cotd(obj)
            out = tvunopebe(@cotd,obj);
        end
        function out = coth(obj)
            out = tvunopebe(@coth,obj);
        end
        function out = csc(obj)
            out = tvunopebe(@csc,obj);
        end
        function out = cscd(obj)
            out = tvunopebe(@cscd,obj);
        end
        function out = csch(obj)
            out = tvunopebe(@csch,obj);
        end
        function out = exp(obj)
            out = tvunopebe(@exp,obj);
        end
        function out = fix(obj)
            out = tvunopebe(@fix,obj);
        end
        function out = fliplr(obj)
            out = tvunopebe(@flip,2,obj);
        end
        function out = flipud(obj)
            out = tvunopebe(@flip,1,obj);
        end
        function out = floor(obj)
            out = tvunopebe(@floor,obj);
        end
        function out = imag(obj)
            out = tvunopebe(@imag,obj);
        end
        function out = isfinite(obj)
            out = tvunopebe(@isfinite,obj);
        end
        function out = isinf(obj)
            out = tvunopebe(@isinf,obj);
        end
        function out = isnan(obj)
            out = tvunopebe(@isnan,obj);
        end
        function out = log(obj)
            out = tvunopebe(@log,obj);
        end
        function out = log10(obj)
            out = tvunopebe(@log10,obj);
        end
        function out = log2(obj)
            out = tvunopebe(@log2,obj);
        end
        function out = not(obj)
            out = tvunopebe(@not,obj);
        end
        function out = real(obj)
            out = tvunopebe(@real,obj);
        end
        function out = round(obj)
            out = tvunopebe(@round,obj);
        end
        function out = sec(obj)
            out = tvunopebe(@sec,obj);
        end
        function out = secd(obj)
            out = tvunopebe(@secd,obj);
        end
        function out = sech(obj)
            out = tvunopebe(@sech,obj);
        end
        function out = sign(obj)
            out = tvunopebe(@sign,obj);
        end
        function out = sin(obj)
            out = tvunopebe(@sin,obj);
        end
        function out = sind(obj)
            out = tvunopebe(@sind,obj);
        end
        function out = sinh(obj)
            out = tvunopebe(@sinh,obj);
        end
        function out = sqrt(obj)
            out = tvunopebe(@sqrt,obj);
        end
        function out = tan(obj)
            out = tvunopebe(@tan,obj);
        end
        function out = tand(obj)
            out = tvunopebe(@tand,obj);
        end
        function out = tanh(obj)
            out = tvunopebe(@tanh,obj);
        end
        function out = uminus(obj)
            out = tvunopebe(@uminus,obj);
        end
        function out = unwrap(varargin)
            out = tvunopebe(@unwrap,varargin{:});
        end
        function out = uplus(obj)
            out = tvunopebe(@uplus,obj);
        end
        
        %--------------Unary Matrix Operations: For-loop over time
        function out = cond(obj,varargin)
            out = tvunop(@cond,obj,varargin{:});
        end
        function out = ctranspose(obj)
            out = tvunop(@ctranspose,obj);
        end
        function out = det(obj)
            out = tvunop(@det,obj);
        end
        function out = diag(obj,varargin)
            out = tvunop(@diag,obj,varargin{:});
        end
        function out = expm(obj)
            out = tvunop(@expm,obj);
        end
        function out = inv(obj)
            out = tvunop(@inv,obj);
        end
        function out = mean(obj,varargin)
            out = tvunop(@mean,obj,varargin{:});
        end
        function out = median(obj,varargin)
            out = tvunop(@median,obj,varargin{:});
        end
        function out = norm(obj,varargin)
            out = tvunop(@norm,obj,varargin{:});
        end
        function out = null(obj,varargin)
            out = tvunop(@null,obj,varargin{:});
        end
        function out = permute(obj,varargin)
            out = tvunop(@permute,obj,varargin{:});
        end
        function out = pinv(obj,varargin)
            out = tvunop(@pinv,obj,varargin{:});
        end
        function out = repmat(obj,varargin)
            out = tvunop(@repmat,obj,varargin{:});
        end
        function out = rank(obj,varargin)
            out = tvunop(@rank,obj,varargin{:});
        end
        function out = reshape(obj,varargin)
            out = tvunop(@reshape,obj,varargin{:});
        end
        function [out,t] = reshapeByTime(obj)
            [m,n] = size(obj);
            [mt,nt] = size(obj.Time);
            out = zeros(max(mt,nt),numel(size(obj.Data(:,:,1))));
            t = obj.Time;
            k = 1;
            for i = 1:m
                for j = 1:n
                    out(:,k) = reshape(obj.Data(i,j,:),max(mt,nt),1);  
                    k = k+1;
                end
            end
        end
        function out = rcond(obj)
            out = tvunop(@rcond,obj);
        end
        function out = squeeze(obj)
            out = tvunop(@squeeze,obj);
        end
        function out = trace(obj)
            out = tvunop(@trace,obj);
        end
        function out = transpose(obj)
            out = tvunop(@transpose,obj);
        end
        
        %--------------Unary Matrix Operations: For-loop + Variable Outputs
        function varargout = chol(obj,varargin)
            varargout = cell(1,nargout);
            [varargout{:}] = tvunop(@chol,obj,varargin{:});
        end
        function varargout = damp(obj)
            varargout = cell(1,nargout);
            [varargout{:}] = tvunop(@damp,obj);
        end
        function varargout = logm(obj)
            varargout = cell(1,nargout);
            [varargout{:}] = tvunop(@logm,obj);
        end
        function varargout = schur(obj,varargin)
            varargout = cell(1,nargout);
            [varargout{:}] = tvunop(@schur,obj,varargin{:});
        end
        function varargout = sqrtm(obj)
            varargout = cell(1,nargout);
            [varargout{:}] = tvunop(@sqrtm,obj);
        end
        
        %--------------Binary Matrix Operations: Single Output
        function out = and(A,B)
            out = tvbinop(@and,A,B);
        end
        function out = atan2(A,B)
            out = tvbinop(@atan2,A,B);
        end
        function out = atan2d(A,B)
            out = tvbinop(@atan2d,A,B);
        end
        function out = eq(A,B)
            out = tvbinop(@eq,A,B);
        end
        function out = ge(A,B)
            out = tvbinop(@ge,A,B);
        end
        function out = gt(A,B)
            out = tvbinop(@gt,A,B);
        end
        function out = ldivide(A,B)
            out = tvbinop(@ldivide,A,B);
        end
        function out = le(A,B)
            out = tvbinop(@le,A,B);
        end
        function out = lft(A,B,varargin)
            out = tvbinop(@lft,A,B,varargin{:});
        end
        function out = lt(A,B)
            out = tvbinop(@lt,A,B);
        end
        function out = minus(A,B)
            out = tvbinop(@minus,A,B);
        end
        function out = mldivide(A,B)
            out = tvbinop(@mldivide,A,B);
        end
        function out = mpower(A,B)
            out = tvbinop(@mpower,A,B);
        end
        function out = mrdivide(A,B)
            out = tvbinop(@mrdivide,A,B);
        end
        function out = mtimes(A,B)
            out = tvbinop(@mtimes,A,B);
        end
        function out = ne(A,B)
            out = tvbinop(@ne,A,B);
        end
        function out = or(A,B)
            out = tvbinop(@or,A,B);
        end
        function out = plus(A,B)
            out = tvbinop(@plus,A,B);
        end
        function out = power(A,B)
            out = tvbinop(@power,A,B);
        end
        function out = rdivide(A,B)
            out = tvbinop(@rdivide,A,B);
        end
        function out = times(A,B)
            out = tvbinop(@times,A,B);
        end
        function out = xor(A,B)
            out = tvbinop(@xor,A,B);
        end
        
        %--------------N-ary Matrix Operations
        function varargout = care(varargin)
            varargout = cell(1,nargout);
            [varargout{:}] = tvnop(@care,varargin{:});
        end
        function varargout = dare(varargin)
            varargout = cell(1,nargout);
            [varargout{:}] = tvnop(@dare,varargin{:});
        end
        function varargout = dlqr(varargin)
            varargout = cell(1,nargout);
            [varargout{:}] = tvnop(@dlqr,varargin{:});
        end
        function  out = dlyap(varargin)
            out = tvnop(@dlyap,varargin{:});
        end
        function  out = dlyapchol(varargin)
            out = tvnop(@dlyapchol,varargin{:});
        end
        function varargout = eig(varargin)
            varargout = cell(1,nargout);
            [varargout{:}] = tvnop(@eig,varargin{:});
        end
        function varargout = lqr(varargin)
            varargout = cell(1,nargout);
            [varargout{:}] = tvnop(@lqr,varargin{:});
        end
        function varargout = lu(varargin)
            varargout = cell(1,nargout);
            [varargout{:}] = tvnop(@lu,varargin{:});
        end
        function  out = lyap(varargin)
            out = tvnop(@lyap,varargin{:});
        end
        function  out = lyapchol(varargin)
            out = tvnop(@lyapchol,varargin{:});
        end
        function varargout = qr(varargin)
            varargout = cell(1,nargout);
            [varargout{:}] = tvnop(@qr,varargin{:});
        end
        function varargout = svd(varargin)
            varargout = cell(1,nargout);
            [varargout{:}] = tvnop(@svd,varargin{:});
        end
        function  out = sylvester(varargin)
            out = tvnop(@sylvester,varargin{:});
        end
        function out = tril(varargin)
            out = tvnop(@tril,varargin{:});
        end
        function out = triu(varargin)
            out = tvnop(@triu,varargin{:});
        end
        
        %--------------Recursive Binary Operations
        function out = blkdiag(varargin)
            out = tvrecop(@blkdiag,varargin{:});
        end
        function out = horzcat(varargin)
            out = tvrecop(@horzcat,varargin{:});
        end
        function out = vertcat(varargin)
            out = tvrecop(@vertcat,varargin{:});
        end
        
        
    end % end of methods
    
    
end % end of classdef




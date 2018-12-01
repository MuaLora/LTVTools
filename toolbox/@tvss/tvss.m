% TVSS Create a time-varying state-space system
%
%   G = tvss(Data,Time) creates a time-varying system. Time is an Nt-by-1
%   vector of increasing time values. If Data is a 1-by-Nt SS array then G
%   is a single-input, single-output time-varying system with value Data(i)
%   at Time(i). If Data is N1-by-N2-by-Nt then G is an N1-by-N2 time-
%   varying system with value Data(:,:,i) at Time(i). In general, if Data
%   is an N1-by-N2-by-...-by-Nad-by-Nt array then G is an N1-by-N2-by-...
%   -by-Nad time-varying SS array with value Data(:,...,:,i) at Time(i).
%
%   G = tvss(Data,Time,InterpolationMethod) specifies the method to 
%   interpolate the state-space matrices between the specified time grid 
%   points.  The method can be 'Linear', 'Flat', 'Nearest', or 'Spline'.
%   Default is InterpolationMethod = 'Linear'.
%
%   G = tvss(Data) generates a constant system.
%
%   G = tvss(A,B,C,D) creates a continuous-time, time-varying system with
%   state matrices defined by TVMATs (A,B,C,D). Additional model properties
%   can be set using name/value pairs.  See help for SS command.
%
%   G = tvss(A,B,C,D,Ts) creates a discrete-time, time-varying system with
%   sample time Ts.
%
%   % EXAMPLE: (CUT/PASTE)
%   % Create 1-state system with time-varying "A" matrix.
%   Time = linspace(0,5,10)';
%   AData = -5+0.1*Time.^2;
%   A = tvmat(AData,Time);
%   B = 1; C = 1; D=0;
%   G = tvss(A,B,C,D);
%
% See also:


classdef (InferiorClasses={?frd, ?ss,?tf,?zpk,?ureal,?ucomplex,?ucomplexm,...
        ?ultidyn,?udyn,?umat,?uss,?ufrd,?tvmat}) tvss
    % Class definition for time-varying system
    
    properties
        Data = [];
        Time = [-inf,inf];
        InterpolationMethod = 'Linear';  % Flat, Nearest, Spline,
        % XXX Add: TimeUnits, Offset data (e.g. as returned by linearize)
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
        % XXX Currently disallows descriptor data.
        function obj = tvss(Data,Time,InterpolationMethod,varargin)
            nin = nargin;
            if nin==0
                % Return default object
                return
            elseif nin==1 
                if isa(Data,'tvss')
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
                % tvss(A,B,C,D,Ts,...) where (A,B,C,D) are TVMATs
                [isConstant,A,B,C,D] = tvswitchyard(Data,Time,...
                    InterpolationMethod,varargin{1});
                varargin = varargin(2:end);               
                
                if isConstant
                   obj.Data = ss(A.Data,B.Data,C.Data,D.Data,varargin{:});
                   obj.isConstant = true;
                else
                   obj.Data = ss(A.Data,B.Data,C.Data,D.Data,varargin{:});
                   obj.Time = A.Time;
                   obj.InterpolationMethod = A.InterpolationMethod;
                end
                
                %error('Invalid syntax for the "tvss" command.');
            end
            obj.Data = ss(obj.Data);
            Nt = numel(obj.Time);
            if ~obj.isConstant
                % XXX-Undocumented: The next code allows Data to be an 
                % Ny-by-Nu-by-Nt-by-1 or an Ny-by-Nu-by-1-by-Nt vector 
                % where Nt is the dim of Time. It reshapes Data to be Ny-
                % -by-Nu-by-1-by-Nt. This ensures that obj.Data always has
                % at least 4 dims with the last dim corresponding to Time. 
                % This simplifies some complexity as Matlab treats an SS
                % array of size Ny-by-Nu-by-Nt-by-1 as having four dims,
                % i.e. it does not drop the trailing "1". As a result, 
                % SS arrays have dim of 2 or >=4 but not 3. 
                if nmodels(obj.Data)==Nt
                    obj.Data = reshape(obj.Data,[1 Nt]);
                end                
            end
                        
            % Use isvalid to perform all error checking
            [pflag,errstr] = isvalid(obj);
            if pflag==0
                error(errstr);
            end
            
            % XXX Squeeze extra singleton dimensions?            
        end
        
        % isvalid
        function [pflag,errstr] = isvalid(obj)
            % ISVALID Determine if TVSS object is valid.
            errstr = [];
            pflag = 1;
            
            % Check Data Class
            Data = obj.Data;
            if ~isa(Data,'ss') 
                pflag = 0;
                errstr = 'Data must be an ss.';
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
            
            % Check Time/Data consistency
            szData = size(obj.Data);
            if obj.isConstant && ~all(Time==[-inf inf])
                pflag = 0;
                errstr = ['Time vector must be [-inf inf] for "Constant" DataType.'];
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
        
        function out = ss(obj)
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
        
                
        function out = ndims(obj)
            if obj.isConstant
                out = ndims( obj.Data );
            else
                out = ndims( obj.Data ) - 1;
                szData = size(obj.Data);
                % SS Arrays have a min dimension = 4
                if out==3 && szData(end-1)==1
                    % Handles case of Ny-by-Nu-by-1-by-NT
                    out=2;
                elseif out==3
                    % Handles case of Ny-by-Nu-by-Na-by-NT
                    out=4;
                end
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
        
        function out = order(obj)
            % XXX Assumes constant state order across time
            out = order(obj.Data(:,:,1));
        end
        
        %--------------Unary Operations: Element-by-Element
        function out = conj(obj)
            out = tvunopebe(@conj,obj);
        end
        function out = uminus(obj)
            out = tvunopebe(@uminus,obj);
        end
        function out = uplus(obj)
            out = tvunopebe(@uplus,obj);
        end
        
        %--------------Unary Matrix Operations: For-loop over time
        function out = ctranspose(obj)
            out = tvunop(@ctranspose,obj);
        end
        function out = inv(obj)
            out = tvunop(@inv,obj);
        end
        function out = norm(obj,varargin)
            out = tvunop(@norm,obj,varargin{:});
        end
        function out = permute(obj,varargin)
            out = tvunop(@permute,obj,varargin{:});
        end
        function out = repmat(obj,varargin)
            out = tvunop(@repmat,obj,varargin{:});
        end
        function out = repsys(obj,varargin)
            out = tvunop(@repsys,obj,varargin{:});
        end
        function out = reshape(obj,varargin)
            out = tvunop(@reshape,obj,varargin{:});
        end
        function out = transpose(obj)
            out = tvunop(@transpose,obj);
        end
        
        %--------------Unary Matrix Operations: For-loop + Variable Outputs
        
        
        %--------------Binary Matrix Operations: Single Output
        function out = feedback(A,B,varargin)
            out = tvbinop(@feedback,A,B,varargin{:});
        end
        function out = ldivide(A,B)
            out = tvbinop(@ldivide,A,B);
        end
        function out = lft(A,B,varargin)
            out = tvbinop(@lft,A,B,varargin{:});
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
        function out = plus(A,B)
            out = tvbinop(@plus,A,B);
        end
        function out = rdivide(A,B)
            out = tvbinop(@rdivide,A,B);
        end
        
        %--------------N-ary Matrix Operations
 
       %--------------Recursive Binary Operations
        function out = append(varargin)
            out = tvrecop(@append,varargin{:});
        end
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
    
    % XXX TODO
    % connect, getLoopTransfer, getIOTransfer
    
    
end % end of classdef




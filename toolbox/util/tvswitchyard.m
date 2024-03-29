function [isConstant,varargout] = tvswitchyard(varargin)
% TVSWITCHYARD Utility function for LTV object interaction.
%
% This function checks all inputs for consistency and lifts objects
% to be TVMAT/TVSS. Specifically, the function:
%  1) Raises empty, double, or logical inputs to be a constant TVMAT/TVSS.  
%  2) Checks that all inputs have a consistent Time & InterpolationMethod. 
%  3) If inputs contain both Constant/Varying objects then the Constant
%     objects are mapped to a Varying object on the correct time vector.
%
% Note: Gridded objects combinations with different time vectors are not
% allowed. Use EVALT to convert to a consistent time vector.

%% Get Time and DataTypes
% Assume each input is either a tvss, tvmat, empty, double, or logical. 
% Empty,  double, or logical are lifted to a Constant tvmat/tvss.
nin = nargin;
varargout = cell(nin,1);
AllTime = cell(nin,1);
AllIM = cell(nin,1);
AllisConstant = false(nin,1);
for i=1:nin
    % Lift empty/double/logical to a TVMAT or TVSS. This will error out
    % if varargin{i} is an input of unexpected class.
    vi = varargin{i};
    if isa(vi,'double') || isa(vi,'logical')
        vi = tvmat(vi);
    elseif isa(vi,'ss') || isa(vi,'tf') || isa(vi,'zpk')
        vi = tvss(vi);
    elseif isa(vi,'umat')
        vi = tvumat(vi);
    elseif isa(vi,'uss')
        vi = tvuss(vi);
    end
     
    varargout{i} = vi;
    AllTime{i} = vi.Time;
    AllIM{i} = vi.InterpolationMethod;
    if vi.isConstant
        AllisConstant(i) = true;
    end
end

%% Check Constant/Varying
% isConstant=true if all inputs are constant
isConstant = all(AllisConstant);

%% Check Time Vectors
if ~isConstant
    DTidx = find(AllisConstant==isConstant);
    Time = AllTime{ DTidx(1) };
    for i=2:numel(DTidx)
        if ~isequal(Time,AllTime{DTidx(i)} )
            error('Operations require a consistent Time vector.');
        end
    end
end

%% Check InterpolationMethod
if ~isConstant
    IM = AllIM{ DTidx(1) };
    for i=2:numel(DTidx)
        if ~isequal(IM,AllIM{DTidx(i)} )
            error('Operations require a consistent InterpolationMethod.');
        end
    end
end

%% Create output
if ~isConstant
    for i=1:nin
        if AllisConstant(i)
            % Lift Constant to TVMAT/TVSS
            Data = varargout{i}.Data;
            nd = ndims(Data)-2;
            Data2 = repmat(Data,[1 1 ones(1,nd) numel(Time)]);
            if isa(Data,'double') || isa(Data,'logical')
                varargout{i} = tvmat(Data2,Time,IM);
            elseif isa(Data,'ss') || isa(Data,'tf') || isa(Data,'zpk')
                varargout{i} = tvss(Data2,Time,IM);
            elseif isa(Data,'umat')
                varargout{i} = tvumat(Data2,Time,IM);
            elseif isa(Data,'uss')
                varargout{i} = tvuss(Data2,Time,IM);
            end
        end
    end
end
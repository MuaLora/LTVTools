classdef (CaseInsensitiveProperties = true, ...
        TruncatedProperties = true) tvodeOptions
    
    % ODE Options set for time varying ODEs
    
    % Examples: 
    %
    % crdeOpt = tvodeOptions;
    % crdeOpt.OdeOptions = odeset('Events',@LOCALevents,'RelTol',1e-5,...
    % 'AbsTol',1e-8);
    % 
    % opt = [];
    
    properties
        % Specify OdeSolver default
        OdeSolver = 'ode45';
        % Option for ODE solver
        OdeOptions = odeset;
    end
    
    methods
        
        %% Specify ODE Solver
        function opt = set.OdeSolver(opt,V)
            % Update this list if new solve becomes available
            availableODESolvers = ...
                {'ode45','ode23','ode113','ode23t','ode23s',...
                'ode23tb','ode15s','ode15i'};
            if any(strcmpi(availableODESolvers, V))
                opt.OdeSolver = lower(V);
            else
                warning('The "OdeSolver" must be set to one of the following:\n%s, %s, %s, %s, %s, %s, %s, %s.\n',...
                    availableODESolvers{:});
            end
        end
        
        %% Specify ODE Options
        function opt = set.OdeOptions(opt,V)
            if isstruct(V) && any(isfield(V,fieldnames(odeset)))
                opt.OdeOptions = V;
            else
                error('The "OdeOptions" must be a structure with fields specified by odeset.');
            end
        end
        
    end
end

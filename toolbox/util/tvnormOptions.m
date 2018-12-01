classdef (CaseInsensitiveProperties = true, ...
        TruncatedProperties = true) tvnormOptions < tvodeOptions
    
    % Options set for TVNORM
    
    properties
        % Display progress of computation [{'off'} | 'on'].
        Display = 'off';       
        % Relative Tolerance for bisection
        RelTol = 1e-3;
        % Absolute Tolerance for bisection
        AbsTol = 1e-4;
        % Initial bounds for bisection
        Bounds = [0, inf];
    end
    
    methods
        function opt = set.Display(opt,V)
            V = ltipack.matchKey(V,{'off','on'});
            if isempty(V)
                error('The "Display" option must be set to ''on'' or ''off''.')
            end
            opt.Display = V;
        end
        
        function opt = set.RelTol(opt,V)
            if isa(V,'double') && isscalar(V) && V >= 0
                opt.RelTol = V;
            else
                error('The "RelTol" option must be a non-negative scalar.')
            end
        end
        
        function opt = set.AbsTol(opt,V)
            if isa(V,'double') && isscalar(V) && V >= 0
                opt.AbsTol = V;
            else
                error('The "AbsTol" option must be a non-negative scalar.')
            end
        end
        
        function opt = set.Bounds(opt,V)
            if isa(V,'double') && numel(V)==2 && V(2) >= V(1) && ...
                    V(1)>=0 && isfinite(V(1))
                opt.Bounds = V;
            else
                error(['The "Bounds" option must be 1-by-2 vector' ...
                    ' with Bounds(2)>=Bounds(1)>=0'])
            end
        end
        
    end
end

classdef Logger < handle
    %LOGGER   Simple logging class with levels and built‑in warning/error
    %
    %   Levels (in increasing severity): DEBUG < INFO < WARN < ERROR
    %
    %   log = Logger();           % create logger
    %   log.setLevel('DEBUG');    % set minimum level to DEBUG
    %
    %   log.debug(  'x = %g', x);                  % prints if level ≤ DEBUG
    %   log.info(   'Loaded %d records', n);       % prints if level ≤ INFO
    %   log.warn(   'MyApp:LowMem', 'Mem = %d MB', m);  % issues warning if level ≤ WARN
    %   log.errorf('MyApp:BadArg', 'Bad arg: %s', name); % throws error
    
    properties (Constant) % read-only
        DEBUG = 0;
        INFO  = 1;
        WARN  = 2;
        ERROR = 3;
    end
    
    properties
        % Current threshold: only messages at or above this will fire
        lvl = Logger.INFO;
    end
    
    methods
        % Constructor
        function obj = Logger(initialLevel)
            %LOGGER Construct an instance, optionally setting its level:
            %  log = Logger()               % defaults to INFO
            %  log = Logger('DEBUG')        % sets threshold to DEBUG
            %  log = Logger(Logger.WARN)    % or pass numeric level
            if nargin>0
                obj.set_level(initialLevel);
            end
        end

        function set_level(obj, lvl)
            %SETLEVEL  Change the threshold by name or numeric value
            if ischar(lvl) || isstring(lvl)
                switch upper(string(lvl))
                    case "DEBUG", obj.lvl = obj.DEBUG;
                    case "INFO",  obj.lvl = obj.INFO;
                    case "WARN",  obj.lvl = obj.WARN;
                    case "ERROR", obj.lvl = obj.ERROR;
                    otherwise
                        error('Logger:BadLevel', 'Unknown log level "%s".', lvl);
                end
            elseif isnumeric(lvl) && isscalar(lvl) && (mod(lvl,1) == 0)
                obj.lvl = lvl;
            else
                error('Logger:BadLevel', 'Invalid level specification.');
            end
        end
        
        function debug(obj, msg, varargin)
            if obj.lvl <= obj.DEBUG
                timestamp = datetime("now");
                fprintf('[%s] DEBUG:\n%s\n', timestamp, sprintf(msg, varargin{:}));
            end
        end
        
        function info(obj, msg, varargin)
            if obj.lvl <= obj.INFO
                timestamp = datetime("now");
                fprintf('Info\n[%s]: %s\n', timestamp, sprintf(msg, varargin{:}));
            end
        end
        
        function warning(obj, id, msg, varargin)
            if obj.lvl <= obj.WARN
                timestamp = datetime("now");
                msg = sprintf('Warn\n[%s]: %s\n', timestamp, sprintf(msg, varargin{:}));
                if nargin >= 3 && ~isempty(id)
                    % NOTE: `builtin()` prevents the object from
                    % recursively calling the its method called warning.
                    % Instead, the builtin function `warning()` is called
                    builtin('warning', id, '%s', msg);
                else
                    builtin('warning', '%s', msg);
                end
            end
        end
        
        function error(~, id, msg, varargin)
            %ERROR  Log a fatal error (always fires, regardless of Level)
            timestamp = datetime("now");
            msg = sprintf('\n[%s] ERROR:\n%s', timestamp, sprintf(msg, varargin{:}));
            if nargin >= 3 && ~isempty(id)
                builtin('error', id, '%s', msg);
            else
                builtin('error', '%s', msg);
            end
        end
    end
end

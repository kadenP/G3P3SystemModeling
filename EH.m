classdef EH < matlab.System & matlab.system.mixin.CustomIcon
    % This MATLAB system block splits a single particle flow inlet into two outlets

    % Public, nontunable properties
    properties (Nontunable) 
        cp_s = 1250                 % (J/kgK) particle specific heat 
    end

    % properties that shouldn't be set by user
    properties (Access = protected) 
      
    end

    % Pre-computed constants
    properties(Access = private)

    end

    methods
        % Constructor
        function obj = TES(varargin)
            % Support name-value pair arguments when constructing object
            setProperties(obj,nargin,varargin{:})
        end
    end
    
    methods(Static, Access = protected)
       %% system block input/output customization
       function icon = getIconImpl(~)
          icon = sprintf('Electric\nHeater'); 
       end
%         function icon = getIconImpl(~)
%             % Define icon for System block
%             icon = matlab.system.display.Icon('png-transparent-simulink-matlab-mathworks-computer-software-logo-coder-miscellaneous-angle-rectangle.png');
%         end
        function [in1name, in2name, in3name] = getInputNamesImpl(~)
          in1name = 'Ts_in';
          in2name = 'mdot_in';
          in3name = 'Tset';
        end
        function [out1name, out2name, out3name] = getOutputNamesImpl(~)
          out1name = 'Ts_out';
          out2name = 'mdot_out';
          out3name = 'Qin';
        end   
       
    end

    methods(Access = protected)
        %% Common functions
        function setupImpl(obj)
            % Perform one-time calculations, such as computing constants 
          
        end
        function resetImpl(obj)
            % Initialize / reset discrete-state properties
            % recomputes static variables. needs to be run if any of the
            % system properties are changed externally.

        end
        function [Ts_out, mdot_out, Qin] = ...
                stepImpl(obj, Ts_in, mdot_in, Tset)
            % Implement algorithm. Calculate y as a function of input u and
            % discrete states.
                        
            % compute outlet flow rates 
            mdot_out = mdot_in;
            if mdot_in > 0 && Ts_in < Tset
                Qin = obj.cp_s*mdot_in*(Tset - Ts_in);
                Ts_out = Tset;
            else
                Ts_out = Ts_in;
                Qin = 0;
            end
        end       
        %% Backup/restore functions
        function s = saveObjectImpl(obj)
            % Set properties in structure s to values in object obj

            % Set public properties and states
            s = saveObjectImpl@matlab.System(obj);

            % Set private and protected properties
            %s.myproperty = obj.myproperty;
        end        
        function loadObjectImpl(obj,s,wasLocked)
            % Set properties in object obj to values in structure s

            % Set private and protected properties
            % obj.myproperty = s.myproperty; 

            % Set public properties and states
            loadObjectImpl@matlab.System(obj,s,wasLocked);
        end
        %% Advanced functions
        function validateInputsImpl(obj,u)
            % Validate inputs to the step method at initialization
        end
        function validatePropertiesImpl(obj)
            % Validate related or interdependent property values
        end
        function processTunedPropertiesImpl(obj)
            % Perform actions when tunable properties change
            % between calls to the System object
        end
        function flag = isInputSizeMutableImpl(obj,index)
            % Return false if input size cannot change
            % between calls to the System object
            flag = false;
        end
        function flag = isInactivePropertyImpl(obj, prop)
            % Return false if property is visible based on object 
            % configuration, for the command line and System block dialog
            flag = false;
        end        
                               
    end
end

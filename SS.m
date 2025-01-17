classdef SS < matlab.System & matlab.system.mixin.CustomIcon
    % This MATLAB system block splits a single particle flow inlet into two outlets

    % Public, nontunable properties
    properties (Nontunable) 
                         
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
          icon = sprintf('Solid\nSplitter'); 
       end
%         function icon = getIconImpl(~)
%             % Define icon for System block
%             icon = matlab.system.display.Icon('png-transparent-simulink-matlab-mathworks-computer-software-logo-coder-miscellaneous-angle-rectangle.png');
%         end
        function [in1name, in2name, in3name] = getInputNamesImpl(~)
          in1name = 'Ts_in';
          in2name = 'mdot_in';
          in3name = 'y1';
        end
        function [out1name, out2name, out3name, out4name] = getOutputNamesImpl(~)
          out1name = 'Ts_out1';
          out2name = 'Ts_out2';
          out3name = 'mdot_out1';
          out4name = 'mdot_out2';
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
        function [Ts_out1, Ts_out2, mdot_out1, mdot_out2] = ...
                stepImpl(obj, Ts_in, mdot_in, y1)
            % Implement algorithm. Calculate y as a function of input u and
            % discrete states.
                        
            % compute outlet flow rates
            Ts_out1 = Ts_in; Ts_out2 = Ts_in;
            mdot_out1 = y1*mdot_in;
            mdot_out2 = (1 - y1)*mdot_in;                                                                                                                                                
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

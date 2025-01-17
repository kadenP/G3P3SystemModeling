classdef HF < matlab.System & matlab.system.mixin.CustomIcon
    % This MATLAB system block simulates the reflected solar radiation from the heliostat field as a function of the incident DNI

    % Public, nontunable properties
    properties (Nontunable) 
        n = 150                     % number of heliostats
        H = 5                       % (m) height of heliostat
        W = 5                       % (m) width of heliostat
        eta_op = 0.6                % opperational efficiency                         
    end

    % properties that shouldn't be set by user
    properties (Access = protected) 
        Aref                        % (m2) area of each heliostat        
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
          icon = sprintf('Falling\nParticle\nReciever'); 
       end
%         function icon = getIconImpl(~)
%             % Define icon for System block
%             icon = matlab.system.display.Icon('png-transparent-simulink-matlab-mathworks-computer-software-logo-coder-miscellaneous-angle-rectangle.png');
%         end
        function [in1name] = getInputNamesImpl(~)
          in1name = 'DNI';
        end
        function [out1name] = getOutputNamesImpl(~)
          out1name = 'Q_rec';
        end   
       
    end

    methods(Access = protected)
        %% Common functions
        function setupImpl(obj)
            % Perform one-time calculations, such as computing constants 
            obj.Aref = obj.H*obj.W;          
        end
        function resetImpl(obj)
            % Initialize / reset discrete-state properties
            % recomputes static variables. needs to be run if any of the
            % system properties are changed externally.

        end
        function [Q_rec] = stepImpl(obj, DNI)
            % Implement algorithm. Calculate y as a function of input u and
            % discrete states.
                        
            % compute reflected solar radiation
            Q_rec = obj.n*obj.Aref*DNI*obj.eta_op;
                                                                                                                                                
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

classdef LSB < matlab.System & matlab.system.mixin.CustomIcon
    % This MATLAB system block simulates transient behavior lumped storage components

    % Public, nontunable properties
    properties (Nontunable)       
        H = 1.5                     % (m) bin height
        D = 0.6                     % (m) bin diameter 
        ms0 = 1                     % (kg) initial mass in bin
        Ts0 = 25                    % (°C) initial temperature of particles
        ht = 10                     % (W/m2K) top surface heat transfer coefficient
        hinf = 10                   % (W/m2K) ambient heat transfer coefficient
        cp_s = 1250                 % (J/kgK) particle specific heat
        rho_s = 3500                % (kg/m3) particle density
        phi_s = 0.6                 % solid volume fraction                                    
    end

    % properties that shouldn't be set by user
    properties (Access = protected)  
        tNow                        % (s) current time
        dt                          % (s) current time step
        As                          % (m2) exposed surface area of bin
        Vbin                        % (m3) volume of hopper 
        xFill                       % particle height fraction
        ms                          % (kg) mass inventory in bin
        
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
          icon = sprintf('Mass\nFlow\nHopper'); 
       end
%         function icon = getIconImpl(~)
%             % Define icon for System block
%             icon = matlab.system.display.Icon('png-transparent-simulink-matlab-mathworks-computer-software-logo-coder-miscellaneous-angle-rectangle.png');
%         end
        function [in1name, in2name, in3name, in4name] = getInputNamesImpl(~)
          in1name = 'Ts_in';
          in2name = 'mdot_s_in';
          in3name = 'mdot_s_out';
          in4name = 't';
        end
        function [out1name, out2name] = getOutputNamesImpl(~)
          out1name = 'Ts_out';
          out2name = 'qloss';
        end  
        function groups = getPropertyGroupsImpl
          group1 = matlab.system.display.SectionGroup( ...
              'Title', 'Geometry Parameters', ...
              'PropertyList', {'H', 'D'});          
          group2 = matlab.system.display.SectionGroup( ...
              'Title', 'Heat Transfer and Material Parameters', ...
              'PropertyList', {'ms0', 'Ts0', 'ht', 'hinf', 'cp_s', 'rho_s', 'phi_s'});                            
          groups = [group1, group2];    
       end
       
    end

    methods(Access = protected)
        %% Common functions
        function setupImpl(obj)
            % Perform one-time calculations, such as computing constants 
            obj.tNow = 0;            
            obj.As = pi*obj.D*obj.H;
            obj.Vbin = pi/4*obj.D^2*obj.H;
        end
        function resetImpl(obj)
            % Initialize / reset discrete-state properties
            % recomputes static variables. needs to be run if any of the
            % system properties are changed externally.

        end
        function [Ts_out, qloss] = ...
                stepImpl(obj, Ts_in, mdot_s_in, mdot_s_out, t)
            % Implement algorithm. Calculate y as a function of input u and
            % discrete states.
            obj.dt = t - obj.tNow;
            
            % step model for new lumped outlet temperature based on inlet
            % and outlet mass flow rates
            
            
                      
            % update time
            obj.tNow = obj.tNow + obj.dt;
                                                                                                                
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
        %% model functions

        
    end
end

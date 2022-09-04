classdef Reciever < matlab.System & matlab.system.mixin.CustomIcon
    % This MATLAB system block can be used to dynamically model the hot and cold thermal storage bins as a function of the particle inlet temperature and the mass flow rate into or out of the bin. It outputs the corresponding particle outlet temperature at each time step.
    %
    % NOTE: When renaming the class name Untitled, the file name
    % and constructor name must be updated to use the class name.
    %
    % This template includes most, but not all, possible properties,
    % attributes, and methods that you can implement for a System object.

    % Public, nontunable properties
    properties (Nontunable) 
        Tinf = 25                   % (�C) ambient temperature
        Ts0 = 25                    % (�C) initial temperature of particles
        cp_s = 1250                 % (J/kgK) particle specific heat
        rho_s = 3500                % (kg/m3) particle density
        phi_s = 0.6                 % solid volume fraction
        sigma = 2.67e-8             % (W/m2K4) stephan-boltzman constant
        epsilon_s = 0.88            % emisivity of falling particles
        alpha_s = 0.88              % absorbtivity of falling particles
        H = 1                       % (m) height of apperature
        W = 1                       % (m) width of apperature
        d = 0.1                     % (m) depth of falling particle curtain       
    end

    % properties that shouldn't be set by user
    properties (Access = protected)  
        tNow                        % (s) current time
        dt                          % (s) current time step 
        mdot                        % (kg/s) FPR flow rate
        Vr                          % (m3) falling curtain volume
        Ar                          % (m2) aperature area
        qRad                        % (W) re-radiaded heat from curtain
        Qabs                        % (W) absorbed radiation
        
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
%        function icon = getIconImpl(~)
%           icon = sprintf('Particle\nTES\nBin'); 
%        end
        function icon = getIconImpl(~)
            % Define icon for System block
            icon = matlab.system.display.Icon('png-transparent-simulink-matlab-mathworks-computer-software-logo-coder-miscellaneous-angle-rectangle.png');
        end
        function [in1name, in2name, in3name, in4name] = getInputNamesImpl(~)
          in1name = 'Ts_in';
          in2name = 'mdot_s_in';
          in3name = 'Qsolar';
          in4name = 't';
        end
        function [out1name, out2name] = getOutputNamesImpl(~)
          out1name = 'Ts_out';
          out2name = 'mdot_s_out';
        end   
        function groups = getPropertyGroupsImpl
          group1 = matlab.system.display.SectionGroup( ...
              'Title', 'Geometry Parameters', ...
              'PropertyList', {'H', 'W', 'd'});          
          group2 = matlab.system.display.SectionGroup( ...
              'Title', 'Heat Transfer and Material Parameters', ...
              'PropertyList', {'cp_s', 'rho_s', 'phi_s', 'epsilon_s', ...
              'alpha_s'});                            
          groups = [group1, group2];    
       end
       
    end

    methods(Access = protected)
        %% Common functions
        function setupImpl(obj)
            % Perform one-time calculations, such as computing constants 
            obj.tNow = 0;
            obj.Vr = obj.H*obj.W*obj.d;
            obj.Ar = obj.H*obj.W;
            

        end
        function resetImpl(obj)
            % Initialize / reset discrete-state properties
            % recomputes static variables. needs to be run if any of the
            % system properties are changed externally.

        end
        function [Ts_out, mdot_s_out] = stepImpl(obj, Ts_in, mdot_s_in, Qsolar)
            % Implement algorithm. Calculate y as a function of input u and
            % discrete states.
            obj.dt = t - obj.tNow;
                        
            % compute variables dependent on inputs
            mdot_s_out = mdot_s_in;
            obj.mdot = mdot_s_in;
            
            % calculate current temperature of particles in apperature
            stepTempSolution(obj, Ts_in, Qsolar);

                      
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
        function stepTempSolution(obj, Ts_in, Qsolar)
            % calculates the next lumped temperature solution for the
            % falling particle temperature outlet
            F_ = @(x) obj.alpha_s*Qsolar - ...
                obj.sigma_s*obj.epsilon_s*obj.Ar*(x^4 - ...
                (obj.T2C(obj.Tinf)^4)) + obj.mdot*obj.cp_s*(x - Ts_in);
            
            
            % implement RK45 method to step to next temperature w F_
            
            
        end
        function T = C2K(~, TC)
           T = TC + 273.15; 
        end
        
        
        
     
    end
end

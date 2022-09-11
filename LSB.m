classdef LSB < matlab.System & matlab.system.mixin.CustomIcon
    % This MATLAB system block simulates transient behavior lumped storage components

    % Public, nontunable properties
    properties (Nontunable)       
        H = 1.5                     % (m) bin height
        D = 0.6                     % (m) bin diameter 
        ms0 = 1                     % (kg) initial mass in bin
        Ts0 = 25                    % (°C) initial temperature of particles
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
        function [in1name, in2name, in3name] = getInputNamesImpl(~)
          in1name = 'Ts_in';
          in2name = 'mdot_s_in';
          in3name = 't';
        end
        function [out1name, out2name, out3name, out4name] = getOutputNamesImpl(~)
          out1name = 'Ts_out';
          out2name = 'mdot_s_out';
          out3name = 'Ts';
          out4name = 'x';
        end   
        function groups = getPropertyGroupsImpl
          group1 = matlab.system.display.SectionGroup( ...
              'Title', 'Geometry Parameters', ...
              'PropertyList', {'n', 'H', 'D'});          
          group2 = matlab.system.display.SectionGroup( ...
              'Title', 'Heat Transfer and Material Parameters', ...
              'PropertyList', {'Ts0', 'cp_s', 'rho_s', 'phi_s'});                            
          groups = [group1, group2];    
       end
       
    end

    methods(Access = protected)
        %% Common functions
        function setupImpl(obj)
            % Perform one-time calculations, such as computing constants 
            obj.tNow = 0;
            obj.x0 = obj.Ts0*ones(obj.n, 1);
            obj.delta = obj.H/(obj.n-2);
            obj.pos = linspace(0, obj.H, obj.n);
            obj.Ac = pi/4*obj.D^2;
            obj.Vbin = obj.Ac*obj.H;
        end
        function resetImpl(obj)
            % Initialize / reset discrete-state properties
            % recomputes static variables. needs to be run if any of the
            % system properties are changed externally.

        end
        function [Ts_out, mdot_s_out, Ts, x_] = ...
                stepImpl(obj, Ts_in, mdot_s_in, t)
            % Implement algorithm. Calculate y as a function of input u and
            % discrete states.
            obj.dt = t - obj.tNow;
            x_ = obj.pos;
                        
            % compute variables dependent on inputs
            mdot_s_out = mdot_s_in;
            obj.v_s = mdot_s_in/(obj.rho_s*obj.phi_s*obj.Ac);     
                  
            % construct linear system matrices and iterate
            buildSystemMatrices(obj);
            Ts = iterateTemps(obj, Ts_in, obj.dt);
            
            % calculate remaining outputs
            Ts_out = Ts(end);
                      
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
        function buildSystemMatrices(obj)
            % uses the current system parameters to construct the linear
            % system matrices to compute the state variables for the next
            % time step
            % particle temperature equations
            Omega1 = obj.v_s/obj.delta;
            Omega2 = -obj.v_s/obj.delta;
            obj.A = spdiags([Omega1*ones(obj.n, 1), ...
                Omega2*ones(obj.n, 1)], [-1, 0], obj.n, obj.n); 
            obj.B = zeros(obj.n, 1); obj.B(1) = Omega1;                      
        end
        function Ts = iterateTemps(obj, Ts_in, t)
            % uses the developed semi-discrete heat exchanger model to
            % compute new temperatures for the particles, sCO2, and metal
            % throughout the heat exchanger.
            % first reformulate as linear system with constant input
            u_ = Ts_in;
            b_ = obj.B*u_;
            Ap = [obj.A, eye(obj.n); zeros(obj.n), zeros(obj.n)];
            xx0 = [obj.x0; b_];               
            xx = expm(t*Ap)*xx0;
            % deconstruct to obtain desired solution
            obj.x = xx(1:obj.n);
            % reset initial conditions
            obj.x0 = obj.x; 
            % assign solution to physical parameters
            Ts = obj.x;          
        end
        
    end
end

classdef BE < matlab.System & matlab.system.mixin.CustomIcon
    % This MATLAB system block simulates the heat loss as particles fall through a certain length of ducting. 

    % Public, nontunable properties
    properties (Nontunable)       
        n = 100                     % number of discretizations in domain
        W = 1                       % (m) elevator case width
        L = 2                       % (m) elevator case length
        Wb = 1/17.25                % (m) bucket case width
        Lb = 2/17.25                % (m) bucket case length
        H = 50                      % (m) elevator case height
        tm = 0.006                  % (m) elevator case thickness
        Ts0 = 20                    % (°C) initial temperature of particles
        Tm0 = 20                    % (°C) initial temperature of metal
        cp_s = 1250                 % (J/kgK) particle specific heat
        rho_s = 3500                % (kg/m3) particle density       
        phi_s = 0.6                 % solid volume fraction
        cp_m = 500                  % (J/kgK) metal specific heat
        rho_m = 8000                % (kg/m3) metal density 
        hinf = 10                   % (W/m2K) ambient convection coefficient  
        wallInsulation = {}         % elevator wall insulation info
        hsw = 180                   % (W/m2K) solid-to-wall convection coefficient
    end

    % properties that shouldn't be set by user
    properties (Access = protected)  
        tNow                        % (s) current time
        dt                          % (s) current time step
        v_s                         % (m/s) bulk solids velocity        
        Acs                         % (m2) open cross-sectional area
        Acm                         % (m2) elevator metal cross-sectional area 
        P1                          % (m) elevator inside perimeter
        P2                          % (m) elevator outside perimeter
        pos                         % (m) position vector
        delta                       % (m) mesh size
        x                           % (°C) temperature states
        x0                          % (°C) initial condition 
        xs                          % (°C) particle temperature states
        xs0                         % (°C) IC for particle states
        xm                          % (°C) metal temperature states
        xm0                         % (°C) IC for metal states
        Ts_in                       % (°C) particle inlet temperature
        A                           % (1/s) 2nx2n state matrix
        B                           % (1/s) 2nx2 input matrix
        As                          % (1/s) nx2n state matrix
        Bs                          % (1/s) 2nx1 input matrix 
        Am                          % (1/s) nx2n state matrix
        Bm                          % (1/s) 2nx1 input matrix 
        Ms                          % (1/s) solid kinematic coefficient
        Msm                         % (1/s) solid-metal kinematic coefficient
        Mminf                       % (1/s) metal-ambient kinematic coefficient
        Uinf                        % (W/m2K) overal ambient h.t. coefficient
        
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
          in3name = 'Tinf';
          in4name = 't';
        end
        function [out1name, out2name, out3name, out4name, out5name, ...
                out6name] = getOutputNamesImpl(~)
          out1name = 'Ts_out';
          out2name = 'mdot_s_out';
          out3name = 'Ts';
          out4name = 'Tm';
          out5name = 'x';
          out6name = 'qloss';
        end   
        function groups = getPropertyGroupsImpl
          group1 = matlab.system.display.SectionGroup( ...
              'Title', 'Geometry Parameters', ...
              'PropertyList', {'n', 'W', 'L', 'H', 'tm'});          
          group2 = matlab.system.display.SectionGroup( ...
              'Title', 'Heat Transfer and Material Parameters', ...
              'PropertyList', {'Ts0', 'Tm0', 'cp_s', 'rho_s', 'phi_s', 'cp_m', ...
                               'rho_m', 'hinf', 'hsw'});                            
          groups = [group1, group2];    
       end
       
    end

    methods(Access = protected)
        %% Common functions
        function setupImpl(obj)
            % Perform one-time calculations, such as computing constants 
            computeUinf(obj);
            obj.tNow = 0;
            obj.xs0 = obj.Ts0*ones(obj.n, 1);
            obj.xm0 = obj.Tm0*ones(obj.n, 1);
            obj.x0 = [obj.xs0; obj.xm0];
            obj.delta = obj.H/(obj.n-2);
            obj.pos = linspace(0, obj.H, obj.n);
            obj.Acs = obj.Wb*obj.Lb;
            obj.Acm = (obj.Wb + obj.tm)*(obj.Lb + obj.tm) - obj.Wb*obj.Lb;
            obj.P1 = 2*obj.Wb + 2*obj.Lb;
            obj.P2 = 2*(obj.W + obj.tm) + 2*(obj.L + obj.tm);
            obj.Ms = obj.hsw*obj.P1/(obj.phi_s*obj.rho_s*obj.cp_s*obj.Acs);
            obj.Msm = obj.hsw*obj.P1/(obj.rho_m*obj.cp_m*obj.Acm);
            obj.Mminf = obj.Uinf*obj.P2/(obj.rho_m*obj.cp_m*obj.Acm);
        end
        function resetImpl(obj)
            % Initialize / reset discrete-state properties
            % recomputes static variables. needs to be run if any of the
            % system properties are changed externally.

        end
        function [Ts_out, mdot_s_out, Ts, Tm, x_, qLoss] = ...
                stepImpl(obj, Ts_in, mdot_s_in, Tinf, t)
            % Implement algorithm. Calculate y as a function of input u and
            % discrete states.
            obj.dt = t - obj.tNow;
            x_ = obj.pos;
            obj.v_s = mdot_s_in/(obj.phi_s*obj.rho_s*obj.Acs);
                        
            % compute variables dependent on inputs
            mdot_s_out = mdot_s_in;    
                  
            % construct linear system matrices and iterate
            buildSystemMatrices(obj);
            [Ts, Tm] = iterateTemps(obj, Ts_in, Tinf, obj.dt);
            
            % calculate remaining outputs
            Ts_out = Ts(end);
            qLoss = mdot_s_in*obj.cp_s*(Ts_in - Ts_out)/1000;
                      
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
            Omega2 = -(obj.v_s/obj.delta + obj.Ms);
            Omega3 = obj.Ms;
            obj.As = spdiags([Omega1*ones(obj.n, 1), ...
                Omega2*ones(obj.n, 1), Omega3*ones(obj.n, 1)], [-1, 0, obj.n], obj.n, 2*obj.n); 
            obj.Bs = zeros(2*obj.n, 1); obj.Bs(1) = Omega1;           
            % metal temperature equations
            Lambda1 = obj.Msm;
            Lambda2 = -(obj.Msm + obj.Mminf);
            obj.Am = spdiags([Lambda1*ones(obj.n, 1), ...
                Lambda2*ones(obj.n, 1)], [0, obj.n], obj.n, 2*obj.n);
            obj.Bm = zeros(2*obj.n, 1); obj.Bm(obj.n+1:end) = obj.Mminf;
            % full state and input matrices
            obj.A = [obj.As; obj.Am];
            obj.B = [obj.Bs, obj.Bm];            
        end
        function [Ts, Tm] = iterateTemps(obj, Ts_in, Tinf, t)
            % uses the developed semi-discrete heat exchanger model to
            % compute new temperatures for the particles, sCO2, and metal
            % throughout the heat exchanger.
            % first reformulate as linear system with constant input
            u_ = [Ts_in; Tinf];
            b_ = obj.B*u_;
            Ap = full([obj.A, eye(2*obj.n); zeros(2*obj.n), zeros(2*obj.n)]);
            xx0 = [obj.x0; b_];               
            xx = expm(t*Ap)*xx0;
            % deconstruct to obtain desired solution
            obj.x = xx(1:2*obj.n);
            obj.xs = obj.x(1:obj.n);
            obj.xm = obj.x(obj.n+1:end);
            % reset initial conditions
            obj.x0 = obj.x;
            obj.xs0 = obj.xs;
            obj.xm0 = obj.xm;  
            % assign solution to physical parameters
            Ts = obj.xs;
            Tm = obj.xm;            
        end  
        function computeUinf(obj)
            % computes the overall heat transfer coefficient with the
            % available insulation information
%             N  = size(obj.wallInsulation{:, 1});
            Rtot = 0;
%             for i = 1:N
%                 x1 = obj.wallInsulation{i, 2}(1);
%                 x2 = obj.wallInsulation{i, 2}(2);
%                 ki = obj.wallInsulation{i, 3};
%                 Rtot = Rtot + (x2 - x1)/ki;                               
%             end
            % add convection
            Rtot = Rtot + 1/obj.hinf;
            obj.Uinf = 1/Rtot;                                               
        end
    end
end

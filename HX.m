classdef HX < matlab.System & matlab.system.mixin.CustomIcon
    % This MATLAB system block simulates transient behavior of the particle-to-sCO2 heat exchanger. The particle, sCO2, and metal temperatures are computed at each time step accross the defined discrete mesh.

    % Public, nontunable properties
    properties (Nontunable)       
        n = 20                      % number of discretizations in domain
        t_s = 0.006                 % (m) solid channel width
        t_CO2 = 0.001               % (m) CO2 channel width
        t_m = 0.003                 % (m) metal thicknesss
        H = 1.5                     % (m) heat exchanger height
        W = 0.6                     % (m) heat exchanger width
        N_plate = 33                % number of parallel plates
        h_s = 450                   % (W/m2K) solid-wall heat transfer coefficient 
        h_CO2 = 3000                % (W/m2K) CO2 heat transfer coefficient      
        Ts0 = 25                    % (°C) initial temperature of particles
        Tco20 = 25                  % (°C) initial temperature of sCO2
        Tm0 = 25                    % (°C) initial temperature of metal
        cp_s = 1250                 % (J/kgK) particle specific heat
        rho_s = 3500                % (kg/m3) particle density
        phi_s = 0.6                 % solid volume fraction        
        cp_CO2 = 1250               % (J/kgK)  specific heat
        rho_CO2 = 110               % (kg/m3) CO2 density        
        cp_m = 500                  % (J/kgK) metal specific heat
        rho_m = 8000                % (kg/m3) metal density                              
    end

    % properties that shouldn't be set by user
    properties (Access = protected)  
        tNow                        % (s) current time
        dt                          % (s) current time step
        vs                          % (m/s) solid velocity       
        vCO2                        % (m/s) CO2 velocity
        u                           % (m/s) velocity inputs
        Q_CO2                       % (W) sCO2 heat addition
        pos                         % (m) position vector
        delta                       % (m) mesh size
        x                           % (°C) temperature states
        x0                          % (°C) initial condition
        xs                          % (°C) particle temperature states
        xCO2                        % (°C) sCO2 temperature states
        xm                          % (°C) metal temperature states
        xs0                         % (°C) particle temperature IC
        xCO20                       % (°C) sCO2 temperature IC
        xm0                         % (°C) metal temperature IC
        w                           % (°C) temperature disturbance inputs
        Ts_in                       % (°C) particle inlet temperature
        Tco2_in                     % (°C) sCO2 inlet temperature
        A                           % (1/s) 3nx3n state matrix
        As                          % (1/s) nx3n particle state matrix eqs
        Aco2                        % (1/s) nx3n sCO2 state matrix eqs
        Am                          % (1/s) nx3n metal state matrix eqs
        B                           % (1/s) 3nx2 input matrix
        Bs                          % (1/s) 3nx1 particle input matrix
        Bco2                        % (1/s) 3nx1 sCO2 input matrix
        kappa_s                     % (1/s) particle kinematic coefficient
        kappa_CO2                   % (1/s) sCO2 kinematic coefficient
        kappa_ms                    % (1/s) metal-s kinematic coefficient
        kappa_mCO2                  % (1/s) metal-sCO2 kinematic coefficient
        xsRef                       % (°C) particle state reference
        xCO2Ref                     % (°C) CO2 state reference
        xmRef                       % (°C) metal state reference
        xRef                        % (°C) state reference
        x0Ref                       % (°C) initial state reference
        vsRef                       % (m/s) solid velocity reference
        vCO2Ref                     % (m/s) CO2 velocity reference
        uRef                        % (m/s) velocity input reference
        Ts_in_Ref                   % (°C) solids inlet temp reference
        Tco2_in_Ref                 % (°C) CO2 inlet temp reference
        fRef                        % (°C/s) state rate reference
        wRef                        % disturbance input reference
        xsPrime                     % (°C) particle state reference
        xCO2Prime                   % (°C) CO2 state reference
        xmPrime                     % (°C) metal state reference
        xPrime                      % (°C) state reference
        x0Prime                     % (°C) state reference
        vsPrime                     % (m/s) solid velocity reference
        vCO2Prime                   % (m/s) CO2 velocity reference
        uPrime                      % (m/s) velocity input reference
        Ts_in_Prime                 % (°C) solids inlet temp reference
        Tco2_in_Prime               % (°C) CO2 inlet temp reference
        fPrime                      % (°C/s) state rate reference
        wPrime                      % disturbance input reference
        Jxs                         % linearized nx3n solid state matrix
        JxCO2                       % linearized nx3n CO2 state matrix
        Jxm                         % linearized nx3n metal state matrix
        Jx                          % linearized 3nx3n state matrix
        Jus                         % linearized 3nx1 solid input matrix
        JuCO2                       % linearized 3nx1 CO2 input matrix
        Ju                          % linearized 3nx2 input matrix
        Jws                         % linearized 3nx1 solid dist. matrix 
        JwCO2                       % linearized 3nx1 CO2 dist. matrix 
        Jw                          % linearized 3nx3 dist. matrix 
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
          icon = sprintf('Particle-to-sCO2\nHeat Exchanger'); 
       end
%         function icon = getIconImpl(~)
%             % Define icon for System block
%             icon = matlab.system.display.Icon('png-transparent-simulink-matlab-mathworks-computer-software-logo-coder-miscellaneous-angle-rectangle.png');
%         end
        function [in1name, in2name, in3name, in4name, in5name] = getInputNamesImpl(~)
          in1name = 'Ts_in';
          in2name = 'Tco2_in';
          in3name = 'mdot_s_in';
          in4name = 'mdot_CO2_in';
          in5name = 't';
        end
        function [out1name, out2name, out3name, out4name, out5name, ...
                  out6name, out7name, out8name, out9name] = getOutputNamesImpl(~)
          out1name = 'Ts_out';
          out2name = 'Tco2_out';
          out3name = 'mdot_s_out';
          out4name = 'mdot_CO2_out';
          out5name = 'Ts';
          out6name = 'Tco2';
          out7name = 'Tm';
          out8name = 'Q_CO2';
          out9name = 'x';
        end   
        function groups = getPropertyGroupsImpl
          group1 = matlab.system.display.SectionGroup( ...
              'Title', 'Geometry Parameters', ...
              'PropertyList', {'n', 't_s', 't_CO2', 't_m', 'H', 'W', ...
                               'N_plate'});          
          group2 = matlab.system.display.SectionGroup( ...
              'Title', 'Heat Transfer and Material Parameters', ...
              'PropertyList', {'h_s', 'h_CO2', 'cp_s', 'rho_s', 'phi_s', ...
                               'cp_CO2', 'rho_CO2', 'cp_m', 'rho_m'});                            
          groups = [group1, group2];    
       end
       
    end

    methods(Access = protected)
        %% Common functions
        function setupImpl(obj)
            % Perform one-time calculations, such as computing constants 
            obj.tNow = 0;
            obj.xs0 = obj.Ts0*ones(obj.n, 1);
            obj.xCO20 = obj.Tco20*ones(obj.n, 1);
            obj.xm0 = obj.Tm0*ones(obj.n, 1);
            obj.x0 = [obj.xs0; obj.xCO20; obj.xm0];
            obj.kappa_s = obj.h_s/(obj.cp_s*obj.rho_s*obj.phi_s*obj.t_s);
            obj.kappa_CO2 = obj.h_CO2/(obj.cp_CO2*obj.rho_CO2*obj.t_CO2);
            obj.kappa_ms = obj.h_s/(obj.cp_m*obj.rho_m*obj.t_m);
            obj.kappa_mCO2 = obj.h_CO2/(obj.cp_m*obj.rho_m*obj.t_m);
            obj.delta = obj.H/(obj.n-2);
            obj.pos = linspace(0, obj.H, obj.n);
            obj.xRef = obj.x0;
            obj.x0Ref = obj.x0;
            obj.xsRef = obj.xs0;
            obj.xCO2Ref = obj.xCO20;
            obj.xmRef = obj.xm0;
            obj.Ts_in_Ref = obj.Ts0;
            obj.Tco2_in_Ref = obj.Tco20;

        end
        function resetImpl(obj)
            % Initialize / reset discrete-state properties
            % recomputes static variables. needs to be run if any of the
            % system properties are changed externally.

        end
        function [Ts_out, Tco2_out, mdot_s_out, mdot_CO2_out, Ts, Tco2, ...
                Tm, Q_CO2, Q_s, x_, Ts_out_lin, Tco2_out_lin] = stepImpl(obj, Ts_in, Tco2_in, mdot_s_in, mdot_CO2_in, t)
            % Implement algorithm. Calculate y as a function of input u and
            % discrete states.
            obj.dt = t - obj.tNow;
            x_ = obj.pos;
                                   
            % compute variables dependent on inputs
            mdot_s_out = mdot_s_in;
            mdot_CO2_out = mdot_CO2_in;
            obj.vs = mdot_s_in/ ...
                      (obj.rho_s*obj.phi_s*obj.t_s*obj.N_plate*obj.W);     
            obj.vCO2 = mdot_CO2_in/ ...
                      (obj.rho_CO2*obj.t_CO2*obj.N_plate*obj.W);
                  
            % compute reference inputs
            if isempty(obj.vsRef), obj.vsRef = obj.vs; end
            if isempty(obj.vCO2Ref), obj.vCO2Ref = obj.vCO2; end
                  
            % construct system matrices and iterate
            buildSystemMatrices(obj);
            [Ts, Tco2, Tm] = iterateTemps(obj, Ts_in, Tco2_in, obj.dt);
            
            % construct linearized system matrices and iterate
            buildLinSystemMatrices(obj);
            [Ts_lin, Tco2_lin, Tm_lin] = iterateLinTemps(obj, Ts_in, Tco2_in, obj.dt);
            Ts_out_lin = Ts_lin(end);
            Tco2_out_lin = Tco2_lin(1);
            
            % calculate remaining outputs
            Ts_out = Ts(end);
            Tco2_out = Tco2(1);
            Q_CO2 = obj.cp_CO2*mdot_CO2_in*(Tco2_out - Tco2_in);
            Q_s = obj.cp_s*mdot_s_in*(Ts_out - Ts_in);
            
            % set new linearization variables
            obj.Ts_in_Ref = Ts_in;
            obj.Tco2_in_Ref = Tco2_in;
                      
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
            Omega1 = obj.vs/obj.delta;
            Omega2 = -(obj.vs/obj.delta + 2*obj.kappa_s);
            Omega3 = 2*obj.kappa_s;
            obj.As = spdiags([Omega1*ones(obj.n, 1), ...
                Omega2*ones(obj.n, 1), Omega3*ones(obj.n, 1)], [-1, 0, 2*obj.n], obj.n, 3*obj.n); 
            obj.Bs = zeros(3*obj.n, 1); obj.Bs(1) = Omega1;
            % sCO2 temperature equations
            Phi1 = -(obj.vCO2/obj.delta + 2*obj.kappa_CO2);
            Phi2 = obj.vCO2/obj.delta;
            Phi3 = 2*obj.kappa_CO2;
            obj.Aco2 = spdiags([Phi1*ones(obj.n, 1), ...
                Phi2*ones(obj.n, 1), Phi3*ones(obj.n, 1)], [obj.n, obj.n+1, 2*obj.n], obj.n, 3*obj.n);
            obj.Aco2(obj.n, 2*obj.n+1) = 0; % hold for sCO2 inlet temp
            obj.Bco2 = zeros(3*obj.n, 1); obj.Bco2(2*obj.n) = Phi2;            
            % metal temperature equations
            Lambda1 = obj.kappa_ms;
            Lambda2 = obj.kappa_mCO2;
            Lambda3 = -(obj.kappa_ms + obj.kappa_mCO2);
            obj.Am = spdiags([Lambda1*ones(obj.n, 1), ...
                Lambda2*ones(obj.n, 1), Lambda3*ones(obj.n, 1)], [0, obj.n, 2*obj.n], obj.n, 3*obj.n);
            % full state and input matrices
            obj.A = [obj.As; obj.Aco2; obj.Am];
            obj.B = [obj.Bs, obj.Bco2];            
        end
        function buildLinSystemMatrices(obj)
            % uses the current system parameters to construct the linear
            % system matrices to compute the state variables for the next
            % time step
            % particle temperature equations
            Omega1Ref = obj.vsRef/obj.delta;
            Omega2Ref = -(obj.vsRef/obj.delta + 2*obj.kappa_s);
            Omega3Ref = 2*obj.kappa_s;
            obj.Jxs = spdiags([Omega1Ref*ones(obj.n, 1), ...
                Omega2Ref*ones(obj.n, 1), Omega3Ref*ones(obj.n, 1)], ...
                [-1, 0, 2*obj.n], obj.n, 3*obj.n); 
            obj.Jus = [0; (obj.xsRef(1:end-1) - obj.xsRef(2:end))/obj.delta; ...
                        zeros(2*obj.n, 1)];
            obj.Jws = zeros(3*obj.n, 1); obj.Jws(1) = Omega1Ref;
            % sCO2 temperature equations
            Phi1Ref = -(obj.vCO2Ref/obj.delta + 2*obj.kappa_CO2);
            Phi2Ref = obj.vCO2Ref/obj.delta;
            Phi3Ref = 2*obj.kappa_CO2;
            obj.JxCO2 = spdiags([Phi1Ref*ones(obj.n, 1), ...
                Phi2Ref*ones(obj.n, 1), Phi3Ref*ones(obj.n, 1)], ...
                [obj.n, obj.n+1, 2*obj.n], obj.n, 3*obj.n);
            obj.JxCO2(obj.n, 2*obj.n+1) = 0; % hold for sCO2 inlet temp
            obj.JuCO2 = [zeros(obj.n, 1); (obj.xCO2Ref(1:end-1) - ...
                       obj.xCO2Ref(2:end))/obj.delta; 0; zeros(obj.n, 1)];
            obj.JwCO2 = zeros(3*obj.n, 1); obj.JwCO2(2*obj.n) = Phi2Ref;            
            % metal temperature equations
            Lambda1Ref = obj.kappa_ms;
            Lambda2Ref = obj.kappa_mCO2;
            Lambda3Ref = -(obj.kappa_ms + obj.kappa_mCO2);
            obj.Jxm = spdiags([Lambda1Ref*ones(obj.n, 1), ...
                Lambda2Ref*ones(obj.n, 1), Lambda3Ref*ones(obj.n, 1)], ...
                [0, obj.n, 2*obj.n], obj.n, 3*obj.n);
            % full state and input matrices
            obj.Jx = [obj.Jxs; obj.JxCO2; obj.Jxm];
            obj.Ju = [obj.Jus, obj.JuCO2];
            obj.Jw = [obj.Jws, obj.JwCO2, eye(3*obj.n)];            
        end
        function [Ts, Tco2, Tm] = iterateTemps(obj, Ts_in, Tco2_in, t)
            % uses the developed semi-discrete heat exchanger model to
            % compute new temperatures for the particles, sCO2, and metal
            % throughout the heat exchanger.
            % first reformulate as linear system with constant input
            w_ = [Ts_in; Tco2_in];
            b_ = obj.B*w_;
            Ap = [obj.A, eye(3*obj.n); zeros(3*obj.n), zeros(3*obj.n)];
            xx0 = [obj.x0; b_];               
            xx = expm(t*Ap)*xx0;
            % deconstruct to obtain desired solution
            obj.x = xx(1:3*obj.n);
            obj.xs = obj.x(1:obj.n);
            obj.xCO2 = obj.x(obj.n+1:2*obj.n);
            obj.xm = obj.x(2*obj.n+1:3*obj.n);
            % reset initial conditions
            obj.x0 = obj.x;
            obj.xs0 = obj.xs;
            obj.xCO20 = obj.xCO2;
            obj.xm0 = obj.xm;  
            % assign solution to physical parameters
            Ts = obj.xs;
            Tco2 = obj.xCO2;
            Tm = obj.xm;            
        end  
        function [Ts, Tco2, Tm] = iterateLinTemps(obj, Ts_in, Tco2_in, t)
            % computes particle, sCO2, and metal temperatures using the
            % linearized model (where velocities are treated as inputs)
            % compute reference ramp function
            obj.fRef = obj.A*obj.xRef + obj.B*[obj.Ts_in_Ref; obj.Tco2_in_Ref];
            % compute differences from reference
            obj.x0Prime = obj.x0 - obj.x0Ref;
            obj.vsPrime = obj.vs - obj.vsRef;
            obj.vCO2Prime = obj.vCO2 - obj.vCO2Ref;
            obj.Ts_in_Prime = Ts_in - obj.Ts_in_Ref;
            obj.Tco2_in_Prime = Tco2_in - obj.Tco2_in_Ref;
            % reformulate as linear system with constant inputs
            u_ = [obj.vsPrime; obj.vCO2Prime];
            w_ = [obj.Ts_in_Prime; obj.Tco2_in_Prime; obj.fRef];
            b_ = obj.Ju*u_ + obj.Jw*w_;
            Ap = [obj.Jx, eye(3*obj.n); zeros(3*obj.n), zeros(3*obj.n)];
            xx0 = [obj.x0Prime; b_];               
            xx = expm(t*Ap)*xx0;
            % deconstruct to obtain desired solution
            obj.xPrime = xx(1:3*obj.n);
            obj.xsPrime = obj.xPrime(1:obj.n);
            obj.xCO2Prime = obj.xPrime(obj.n+1:2*obj.n);
            obj.xmPrime = obj.xPrime(2*obj.n+1:3*obj.n);
            obj.x = obj.xPrime + obj.xRef;
            obj.xs = obj.x(1:obj.n);
            obj.xCO2 = obj.x(obj.n+1:2*obj.n);
            obj.xm = obj.x(2*obj.n+1:3*obj.n);
            % reset linearization parameters
            obj.xRef = obj.x0;
            obj.xsRef = obj.xs0;
            obj.xCO2Ref = obj.xCO20;
            obj.xmRef = obj.xm0;
            % reset initial conditions
            obj.x0 = obj.x;
            obj.xs0 = obj.xs;
            obj.xCO20 = obj.xCO2;
            obj.xm0 = obj.xm;  
            % assign solution to physical parameters
            Ts = obj.xs;
            Tco2 = obj.xCO2;
            Tm = obj.xm;
            
            
            
        end
                  
    end
end

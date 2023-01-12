classdef FPR < matlab.System
    % This MATLAB system block simulates the heating of particles falling through the falling particle reciever. It is approximated as a lumped system with the solar flux coming directly from the heliostat field

    % Public, nontunable properties
    properties (Nontunable) 
        Ts0 = 25                    % (°C) initial temperature of particles        
        hInf = 5                    % (W/m2K) ambient heat transfer coefficient
        cp_s = 1250                 % (J/kgK) particle specific heat
        rho_s = 3500                % (kg/m3) particle density
        phi_s = 0.6                 % solid volume fraction
        sigma = 5.67e-8             % (W/m2K4) stephan-boltzman constant
        epsilon_s = 0.88            % emisivity of falling particles
        alpha_s = 0.92              % absorbtivity of falling particles
        H = 1.2                     % (m) height of apperature
        W = 1.2                     % (m) width of apperature
        d = 0.0417                  % (m) depth of falling particle curtain 
        TinfRef = 0                 % (°C) reference ambient temperature
        TinRef = 600                % (°C) reference inlet temperature
        qsRef = 0                   % (W/m2) reference flux (0)        
        dtLin = 0.5                 % (s) maximum time step for lin model 
        tauLag = 10                 % (s) lag filter time constant
        tauLead = 1                 % (s) lead filter time constant
    end

    % properties that shouldn't be set by user
    properties (Access = protected) 
        Tinf                        % (°C) ambient temperature (input)
        Tset                        % (°C) setpoint temperature
        x0                          % (°C) IC for curtain temp
        tNow                        % (s) current time
        dt                          % (s) current time step 
        mdot                        % (kg/s) FPR flow rate
        Vr                          % (m3) falling curtain volume
        Ar                          % (m2) aperature area
        qRad                        % (W) re-radiaded heat from curtain
        Qabs                        % (W) absorbed radiation
        kappaSol                    % (m2-K/J) solar input coefficient
        kappaRad                    % (1/s-K3) re-radiation coefficient
        kappaConv                   % (1/s) convective loss coefficient
        kappaAdv                    % (1/kg) advective loss coefficient
        beta                        % (1/s) linearized state coefficient
        gamma                       % (1/s) linearized ambient coefficient
        theta                       % (1/s) linearized inlet coefficient
        psi                         % (m2-K/J) linearized flux coefficient
        zeta                        % (K/kg) linearized flow coefficient
        TsPrime                     % (K) linearized state
        TinfPrime                   % (K) linearized ambient temperature
        TinPrime                    % (K) linearized inlet temperature
        qsPrime                     % (W/m2) linearized flux
        mdotPrime                   % (kg/s) linearized mass flow rate
        TsRef                       % (K) reference state  
        F_Ref                       % (°C/s) reference rate
        mdotRef                     % (kg/s) reference mass flow rate
        A                           % linearized state matrix
        Bu                          % linearized input matrix
        Bw                          % linearized disturbance vector
        Cz                          % measurment state matrix
        Dzw                         % measurment disturbance matrix
        w                           % linearized disturbance inputs
        Ky                          % feedback controller
        Kuw                         % feedforward controller
        ks                          % scaling vector for controller
                
    end

    % Pre-computed constants
    properties(Access = private)

    end

    methods
        % Constructor
        function obj = TES(varargin)
            % Support name-value pair arguments when constructing object
            setProperties(obj, nargin, varargin{:})
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
        function [in1name, in2name, in3name, in4name, in5name] = getInputNamesImpl(~)
          in1name = 'Ts_in';
          in2name = 'Tinf';
          in3name = 'Qsolar';
          in4name = 't';
          in5name = 'Tset';
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
              'PropertyList', {'Ts0', 'hInf', 'cp_s', 'rho_s', 'phi_s', 'epsilon_s', ...
              'alpha_s'}); 
          group3 = matlab.system.display.SectionGroup( ...
              'Title', 'Control Parameters', ...
              'PropertyList', {'tauLag', 'tauLead'});
          groups = [group1, group2, group3];    
       end
       
    end

    methods(Access = protected)
        %% Common functions
        function setupImpl(obj)
            % Perform one-time calculations, such as computing constants 
            obj.tNow = 0;
            obj.x0 = obj.Ts0;
            obj.mdotRef = 1;
            obj.TsRef = obj.Ts0;
            obj.mdot = 1;
            obj.Vr = obj.H*obj.W*obj.d;
            obj.Ar = obj.H*obj.W;
            obj.kappaSol = obj.alpha_s/(obj.cp_s*obj.rho_s*obj.phi_s*obj.d);
            obj.kappaRad = 1*obj.epsilon_s*obj.sigma/...
                                    (obj.cp_s*obj.rho_s*obj.phi_s*obj.d);
            obj.kappaConv = 2*obj.hInf/(obj.cp_s*obj.rho_s*obj.phi_s*obj.d);
            obj.kappaAdv = 1/(obj.rho_s*obj.phi_s*obj.Vr);
            obj.ks = [1, 1, obj.hInf, 1];
            

        end
        function resetImpl(obj)
            % Initialize / reset discrete-state properties
            % recomputes static variables. needs to be run if any of the
            % system properties are changed externally.

        end
        function [Ts_out, mdot_s_out] = stepImpl(obj, Ts_in, Tinf, Qsolar, t, Tset)
            % Implement algorithm. Calculate y as a function of input u and
            % discrete states.
            obj.dt = t - obj.tNow;
                        
            % compute variables dependent on inputs
            obj.Tinf = Tinf;
            obj.Tset = Tset;            
            
            % generate mass flow control law with linear model
            if Qsolar > 0
                obj.mdot = setMassFlowRate(obj, Ts_in, Qsolar/obj.Ar, obj.dt);
            else
                obj.mdot = 0;
            end
            mdot_s_out = obj.mdot;
            
            % compute temperature at next step
            Ts_out = stepTempSolution(obj, Ts_in, Qsolar/obj.Ar);
                                              
            % reset the reference state
            obj.TsRef = Ts_out;
            obj.mdotRef = obj.mdot;
                                      
            % update time and initial condition
            obj.tNow = obj.tNow + obj.dt;
            obj.x0 = Ts_out;
                                                                                                                
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
        function Ts = stepTempSolution(obj, Ts_in, qsolar)
            % calculates the next lumped temperature solution for the
            % falling particle temperature outlet
            F_ = @(x) obj.kappaSol*qsolar - ...
                obj.kappaRad*(x^4 - (obj.C2K(obj.Tinf)^4)) - ...
                obj.kappaConv*(x - (obj.C2K(obj.Tinf))) - ...
                obj.kappaAdv*obj.mdot*(x - obj.C2K(Ts_in));
                       
            % implement second order midpoint method to step to next temperature w F_
            a = 0.5;
            TsK = obj.C2K(obj.x0) + ...
                obj.dt*((1 - 1/(2*a))*F_(obj.C2K(obj.x0)) + ...
                1/(2*a)*F_(obj.C2K(obj.x0) + a*obj.dt*F_(obj.C2K(obj.x0))));
            Ts = TsK - 273.15;
                        
        end
        function buildLinTempSystem(obj, Ts_in, qsolar)
            % sets Jacobians and system matrices for the linearized FPR model            
            obj.F_Ref = obj.kappaSol*obj.qsRef - ...
                    obj.kappaRad*((obj.C2K(obj.TsRef))^4 - (obj.C2K(obj.TinfRef)^4)) - ...
                    obj.kappaConv*(obj.TsRef - obj.TinfRef) - ...
                    obj.kappaAdv*obj.mdotRef*(obj.TsRef - obj.TinRef);

            % set Jacobian variables
            obj.gamma = 4*obj.kappaRad*(obj.C2K(obj.TinfRef))^3 ...
                        + obj.kappaConv;
            obj.theta = obj.kappaAdv*obj.mdotRef;
            obj.psi = obj.kappaSol;
            obj.beta = -4*obj.kappaRad*(obj.C2K(obj.TsRef))^3 - obj.kappaConv ...
                        - obj.kappaAdv*obj.mdotRef;
            obj.zeta = obj.kappaAdv*(obj.TinRef - obj.TsRef);

            % set linearized state-space variables
            obj.TsPrime = obj.x0 - obj.TsRef;
            obj.TinfPrime = obj.Tinf - obj.TinfRef;
            obj.TinPrime = Ts_in - obj.TinRef;
            obj.qsPrime = qsolar - obj.qsRef;
            obj.mdotPrime = obj.mdot - obj.mdotRef;

            % set linear system matrices
            obj.A = obj.beta;
            obj.Bu = obj.zeta;
            obj.Bw = [obj.gamma, obj.theta, obj.psi, 1];
            obj.Cz = [1; 0; 0; 0];
            obj.Dzw = [0, 0, 0, 0; 1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0];
          
        end
        function Ts = stepLinTempSolution(obj, Ts_in, qsolar)
            % calculates the next lumped temperature solution with the
            % linearized system model
            
            % set Jacobian variables;
            obj.gamma = 4*obj.kappaRad*(obj.C2K(obj.TinfRef))^3 ...
                        + obj.kappaConv;
            obj.theta = obj.kappaAdv*obj.mdotRef;
            obj.psi = obj.kappaSol;
            
            % reduce time step if too large
            if obj.dt > obj.dtLin
                t_ = linspace(0, obj.dt, ceil(obj.dt/obj.dtLin));
                dt_ = t_(2);
            else
                t_ = obj.dt; 
                dt_ = obj.dt;
            end
            x0_ = obj.x0;
            
            for i = 1:length(t_)           
                obj.F_Ref = obj.kappaSol*obj.qsRef - ...
                    obj.kappaRad*((obj.C2K(obj.TsRef))^4 - (obj.C2K(obj.TinfRef)^4)) - ...
                    obj.kappaConv*(obj.TsRef - obj.TinfRef) - ...
                    obj.kappaAdv*obj.mdotRef*(obj.TsRef - obj.TinRef);

                % set Jacobian variables
                obj.beta = -4*obj.kappaRad*(obj.C2K(obj.TsRef))^3 - obj.kappaConv ...
                            - obj.kappaAdv*obj.mdotRef;
                obj.zeta = obj.kappaAdv*(obj.TinRef - obj.TsRef);

                % set state-space variables
                obj.TsPrime = x0_ - obj.TsRef;
                obj.TinfPrime = obj.Tinf - obj.TinfRef;
                obj.TinPrime = Ts_in - obj.TinRef;
                obj.qsPrime = qsolar - obj.qsRef;
                obj.mdotPrime = obj.mdot - obj.mdotRef;
                A_ = obj.beta;
                B = [obj.zeta, obj.gamma, obj.theta, obj.psi, 1];
                u = [obj.mdotPrime; obj.TinfPrime; obj.TinPrime; ...
                    obj.qsPrime; obj.F_Ref];
                b_ = B*u;

                % step linear model solution
                Ap = [A_, eye(size(A_)); zeros(size(A_)), zeros(size(A_))];
                xx0 = [obj.TsPrime; b_];               
                xx = expm(dt_.*Ap)*xx0;
                Ts = xx(1) + obj.TsRef;
                
                % reset interim variables
                obj.TsRef = x0_;
                x0_ = Ts;               
            end                                                         
        end  
        function [mdot_] = setMassFlowRate(obj, Ts_in, qsolar, tg)
            % uses feedback control to set the mass flow rate according to
            % the current temperature error
            
            % set linearized system parameters
            buildLinTempSystem(obj, Ts_in, qsolar)                                         

            % set feedback and feedforward controller
            if abs(obj.zeta) > obj.kappaAdv*10
                % dynamic feedforward controller
                obj.Kuw = -obj.Bw./obj.zeta - ...
                        pinv(obj.Cz)*obj.Dzw.*exp(-obj.beta*tg)./obj.zeta;
                mdot_ff = obj.mdotRef + obj.Kuw*[obj.TinfPrime; obj.TinPrime; ...
                    obj.qsPrime; obj.F_Ref];
                
                % feedback signal with lead-lag compensation               
                Fll = (1 - exp(-(tg)/obj.tauLag))/ ...
                      (1 - exp(-(tg)/obj.tauLead));
                obj.Ky = Fll*obj.ks/obj.zeta*[obj.gamma; obj.theta; obj.psi; 1];
                mdot_fb = obj.Ky*(obj.Tset - obj.x0); 
                
                % feedback and feedforward signals combined
                obj.mdot = mdot_ff + mdot_fb;
            end            
            
            % limiting conditions
            if obj.mdot <= 0, obj.mdot = 0; end
            if obj.mdot >= 10, obj.mdot = 10; end                      
                                                       
            % set function outputs
            mdot_ = obj.mdot;            
                       
        end
        function T = C2K(~, TC)
            T = TC + 273.15; 
        end
           
    end
end

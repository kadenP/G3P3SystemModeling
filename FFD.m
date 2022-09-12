classdef FFD < matlab.System & matlab.system.mixin.CustomIcon
    % This MATLAB system block simulates the heat loss as particles fall through a certain length of ducting. 

    % Public, nontunable properties
    properties (Nontunable)       
        n = 100                     % number of discretizations in domain
        D = 0.25                    % (m) duct inside diameter
        L = 1                       % (m) duct length
        td = 0.003                  % (m) duct thickness
        Ts0 = 25                    % (°C) initial temperature of particles
        Tm0 = 25                    % (°C) initial temperature of metal
        cp_s = 1250                 % (J/kgK) particle specific heat
        rho_s = 3500                % (kg/m3) particle density       
        phi_s = 0.6                 % solid volume fraction
        cp_m = 500                  % (J/kgK) metal specific heat
        rho_m = 8000                % (kg/m3) metal density 
        hinf = 10                   % (W/m2K) ambient convection coefficient    
        hsw = 100                   % (W/m2K) solid-to-wall convection coefficient
    end

    % properties that shouldn't be set by user
    properties (Access = protected)  
        tNow                        % (s) current time
        dt                          % (s) current time step
        v_s                         % (m/s) bulk solids velocity        
        Acs                         % (m2) open cross-sectional area
        Acm                         % (m2) duct cross-sectional area 
        P1                          % (m) duct inside perimeter
        P2                          % (m) duct outside perimeter
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
        function [out1name, out2name, out3name, out4name, out5name] = ...
                getOutputNamesImpl(~)
          out1name = 'Ts_out';
          out2name = 'mdot_s_out';
          out3name = 'Ts';
          out4name = 'Tm';
          out5name = 'x';
        end   
        function groups = getPropertyGroupsImpl
          group1 = matlab.system.display.SectionGroup( ...
              'Title', 'Geometry Parameters', ...
              'PropertyList', {'n', 'D', 'L', 'td'});          
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
            obj.tNow = 0;
            obj.xs0 = obj.Ts0*ones(obj.n, 1);
            obj.xm0 = obj.Tm0*ones(obj.n, 1);
            obj.x0 = [obj.xs0; obj.xm0];
            obj.delta = obj.L/(obj.n-2);
            obj.pos = linspace(0, obj.L, obj.n);
            obj.Acs = pi/4*obj.D^2;
            obj.Acm = pi/4*((obj.D + obj.td)^2 - obj.D^2);
            obj.P1 = pi*obj.D;
            obj.P2 = pi*(obj.D + obj.td);
            obj.Ms = obj.hsw*obj.P1/(obj.phi_s*obj.rho_s*obj.cp_s*obj.Acs);
            obj.Msm = obj.hsw*obj.P1/(obj.rho_m*obj.cp_m*obj.Acm);
            obj.Mminf = obj.hinf*obj.P2/(obj.rho_m*obj.cp_m*obj.Acm);
        end
        function resetImpl(obj)
            % Initialize / reset discrete-state properties
            % recomputes static variables. needs to be run if any of the
            % system properties are changed externally.

        end
        function [Ts_out, mdot_s_out, Ts, Tm, x_] = ...
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
            Ap = [obj.A, eye(2*obj.n); zeros(2*obj.n), zeros(2*obj.n)];
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
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % top boundary temperature distribution for discharge mode
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function thetaT_ = iterateThetaT(obj, thetaS_, r_, topIC, zhat_, rtop_, zcenter_)
            % iterates the temperature for discretized energy equation
            n = length(zhat_); m = length(rtop_);
            % compute velocity distribution if not populated already
            ubar_ = computeUbar(obj, rtop_);
            % compute matrix exponential if empty
            [GT_, expTop_, LambdaT2_] = computeExpTop(obj, ubar_, zhat_, rtop_);
            % boundary condition at z=0 
            topIC(end, :) = obj.S2T(thetaS_, r_, rtop_);            
            % boundary condition at r = b
            topIC(:, end) = thetaS_(end, end);
            % set initial condition           
            topIC = topIC';
            topIC = GT_*topIC(:) + LambdaT2_*obj.thetaA;
            % calculate new temperature states with matrix exponential                    
            thetaT_ = expTop_*topIC(:);
            % set new temperatures for top region
            thetaT_ = reshape(thetaT_, m, n)';
            % fix corner temperatures
            thetaT_(1, 1) = thetaT_(1, 2);
            thetaT_(end, 1) = thetaT_(end, 2);
            thetaT_(1, end) = thetaT_(2, end);
            thetaT_(end, end) = thetaT_(end-1, end);            
            % compute upper boundary temperature for center channel
            computeThetaChat(obj, thetaT_, ubar_, zcenter_);
        end      
        function [GT_, expTop_, LambdaT2_] = computeExpTop(obj, ubar_, zhat_, rtop_)
            % computes the static top exponential
            n = length(zhat_); m = length(rtop_); nm = n*m;
            obj.drtop = rtop_(2) - rtop_(1);
            obj.dzhat = zhat_(2) - zhat_(1);
            % top boundary meshed domain state matrix elements 
            OmegaT1_ = -2/obj.drtop^2 - 2/obj.dzhat^2;
            Ls = diag(ones(length(obj.ubar)-1, 1), -1);
            OmegaT2_ = 1/obj.drtop^2 + (1./repmat(rtop_', n, 1) ... 
                   - obj.Pe*repmat(Ls*ubar_', n, 1))/(2*obj.drtop);
            OmegaT3_ = 1/obj.dzhat^2;
            Us = diag(ones(length(ubar_)-1, 1), 1);
            OmegaT4_ = 1/obj.drtop^2 - (1./repmat(rtop_', n, 1) ... 
                   - obj.Pe*repmat(Us*ubar_', n, 1))/(2*obj.drtop);
            OmegaT5_ = 1/obj.dzhat^2;
            AT_ = spdiags([OmegaT1_*ones(nm, 1), OmegaT2_, ...
                        OmegaT3_*ones(nm, 1), OmegaT4_, ...
                        OmegaT5_*ones(nm, 1)], [0, 1, m, -1, -m], ...
                        nm, nm);                                         
            % eliminate boundary elements from AT
            AT_(1:m, :) = 0;                 % z = h
            AT_(end-m+1:end, :) = 0;         % z = 0
            AT_(m:m:end, :) = 0;             % r = b
            AT_(1:m:end-m+1, :) = 0;         % r = a
            % Boundary state matrix
            GT_ = eye(nm);
            % boundary condition at z=h
            AT_(1:m, :) = AT_(m+1:2*m, :)*1/(1 + obj.Bi5*obj.dzhat); 
            LambdaT2_ = [obj.Bi5*obj.dzhat/(1 + obj.Bi5*obj.dzhat) ...
                              *ones(m-1, 1); zeros(nm-m+1, 1)];
            GT_((m-1)*nm:nm+1:2*m*nm) = 1/(1 + obj.Bi5*obj.dzhat);
            GT_(1:nm+1:m*(nm+1)) = 0;
            % boundary condition at r = a
            AT_(1:m:end-m+1, :) = AT_(2:m:end-m+2, :);
            % matrix exponential for single time step
            expTop_ = expm(AT_*obj.df);
        end
        function computeThetaChat(obj, thetaT_, ubar_, zcenter_)
            % computes the homogeneous temperature for the top boundary in
            % the center channel based on an energy balance on the junction
            % that connects the top and center boundaries
            wbar_ = computeWbar(obj, zcenter_);
            fz = simpsonIntegrator(obj, obj.zhat);
            obj.thetaChat = 2/obj.a0*ubar_(1)/wbar_(1)*fz*thetaT_(:, 1);           
        end
        function recordTopLoss(obj)
            % records dT from outer radius to inner radius in the top
            % channel at the current time step
            if isempty(obj.topLoss); obj.topLoss = cell(3, 1); end
            obj.topLoss{1} = [obj.topLoss{1}, obj.FoNow];
            obj.topLoss{2} = [obj.topLoss{2}, ...
                                  obj.thetaT(:, 2) - obj.thetaT(:, end-1)];
            obj.topLoss{3} = [obj.topLoss{3}, ...
                              obj.theta2T(obj.thetaT(:, 2)) ...
                            - obj.theta2T(obj.thetaT(:, end-1))];
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % center channel velocity distribution for discharge mode
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        function wbar_ = computeWbar(obj, zcenter_)
            % compute bulk velocity distribution in center channel
            wbar_ = obj.Qc/(pi*(obj.a0)^2)*ones(size(zcenter_));             
        end        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % center boundary temperature distribution for discharge mode
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function thetaC_ = iterateThetaC(obj, thetaS_, z_, centerIC, zcenter_, rhat_)
            % iterates the temperature for discretized energy equation
            n = length(zcenter_); m = length(rhat_); nm = n*m; 
            obj.dzc = zcenter_(2) - zcenter_(1);
            obj.drhat = rhat_(2) - rhat_(1);
            % compute velocity distribution if not populated already
            wbar_ = computeWbar(obj, zcenter_);
            % top boundary meshed domain state matrix elements 
            OmegaC1_ = -2/obj.dzc^2 - 2/obj.drhat^2;
            OmegaC2_ = 1/obj.drhat^2 ...
                        + 1/(2*repmat(rhat_', n, 1)*obj.drhat);
            OmegaC3_ = 1/obj.dzc^2 + 0.1*obj.Pe ...
                   *repmat(wbar_', m, 1)/(2*obj.dzc);
            OmegaC4_ = 1/obj.drhat^2 ...
                        - 1/(2*repmat(rhat_', n, 1)*obj.drhat);
            OmegaC5_ = 1/obj.dzc^2 - 0.1*obj.Pe ...
                   *repmat(wbar_', m, 1)/(2*obj.dzc);
            AC_ = spdiags([OmegaC1_*ones(nm, 1), OmegaC2_', ...
                              OmegaC3_, OmegaC4_', OmegaC5_], ...
                              [0, 1, m, -1, -m], nm, nm);
            % eliminate boundary elements from AC
            AC_(1:m, :) = 0;                 % z = 1
            AC_(end-m+1:end, :) = 0;         % z = 0
            AC_(m:m:end, :) = 0;             % r = a
            AC_(1:m:end-m+1, :) = 0;         % r = 0            
            % set boundary temperatures
            [thetaC_, AC_] = setCenterChannelBoundaries(obj, AC_, thetaS_, z_, centerIC, zcenter_, rhat_);
            % compute new temperature states with matrix exponential
            thetaC_ = thetaC_';
            thetaC_ = expm(AC_*obj.df)*thetaC_(:); 
            % set new temperatures for top region
            thetaC_ = reshape(thetaC_, m, n)';
        end           
        function [thetaC_, AC_] = setCenterChannelBoundaries(obj, AC_, thetaS_, z_, thetaC_, zcenter_, rhat_)
            % sets center boundary conditions
            m = length(rhat_);
            % boundary condition at r = 0  
            AC_(1:m:end-m+1, :) = AC_(2:m:end-m+2, :);
            thetaC_(:, 1) = thetaC_(:, 2);
            % boundary condition at z=0
            AC_(end-m+1:end, :) = AC_(end-2*m+1:end-m, :);
            thetaC_(end, 2:end) = thetaC_(end-1, 2:end);
            % boundary condition at r = a
            thetaC_(:, end) = obj.S2C(thetaS_, z_, zcenter_); 
            % boundary condition at z=1
            thetaC_(1, :) = obj.thetaChat;
        end                      
    end
end

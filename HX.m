classdef HX < matlab.System & matlab.system.mixin.CustomIcon
    % This MATLAB system block can be used to dynamically model the hot and cold thermal storage bins as a function of the particle inlet temperature and the mass flow rate into or out of the bin. It outputs the corresponding particle outlet temperature at each time step.
    %
    % NOTE: When renaming the class name Untitled, the file name
    % and constructor name must be updated to use the class name.
    %
    % This template includes most, but not all, possible properties,
    % attributes, and methods that you can implement for a System object.

    % Public, nontunable properties
    properties (Nontunable)       
        n = 100                     % number of discretizations in domain
        hc_s = 0.006                % (m) solid channel width
        hc_CO2 = 0.001              % (m) CO2 channel width
        t_m = 0.003                 % (m) metal thicknesss
        H = 1                       % (m) heat exchanger height
        W = 0.5                     % (m) heat exchanger width
        N_plate = 25                % number of parallel plates
        h_conv_CO2 = 1000           % (W/m2K) CO2 heat transfer coefficient
        h_conv_sw = 180             % (W/m2K) solid-wall heat transfer coefficient        
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
        v_s                         % (m/s) solid velocity       
        v_CO2                       % (m/s) CO2 velocity
        Q_CO2                       % (W) sCO2 heat addition
        pos                         % (m) position vector
        delta                       % (m) mesh size
        x                           % (°C) temperature states
        x0                          % (°C) initial condition
        xs                          % (°C) particle temperature states
        xCO2                        % (°C) sCO2 temperature states
        xm                          % (°C) metal temperature states
        u                           % (°C) temperature inputs
        Ts_in                       % (°C) particle inlet temperature
        Tco2_in                     % (°C) particle inlet temperature
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
          in2name = 'Tco2_in';
          in3name = 'mdot_s_in';
          in4name = 'mdot_CO2_in';
        end
        function [out1name, out2name, out3name, out4name, out5name, ...
                  out6name, out7name, out8name] = getOutputNamesImpl(~)
          out1name = 'Tsout';
          out2name = 'Tco2_out';
          out3name = 'mdot_s_out';
          out4name = 'mdot_CO2_out';
          out5name = 'Ts';
          out6name = 'Tco2';
          out7name = 'Tm';
          out8name = 'Q_CO2';
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

        end
        function resetImpl(obj)
            % Initialize / reset discrete-state properties
            % recomputes static variables. needs to be run if any of the
            % system properties are changed externally.

        end
        function [y1, y2, y3, y4, y5, y6, y7, y8] = ...
                stepImpl(obj, Ts_in, Tco2_in, mdot_s_in, mdot_CO2_in)
            % Implement algorithm. Calculate y as a function of input u and
            % discrete states.
           
            
            % update initial condition
            obj.FoNow = obj.FoNow + obj.df;
            obj.IC(:) = NaN; obj.rIC(:) = NaN; obj.zIC(:) = NaN;
            [n_, m_] = size(IC_);
            obj.IC(1:n_, 1:m_) = IC_;
            obj.zIC(1:n_) = z_;
            obj.rIC(1:m_) = r_;
                                                                                                                
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

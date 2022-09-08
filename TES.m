classdef TES < matlab.System & matlab.system.mixin.CustomIcon
    % This MATLAB system block can be used to dynamically model the hot and cold thermal storage bins as a function of the particle inlet temperature and the mass flow rate into or out of the bin. It outputs the corresponding particle outlet temperature at each time step.
    %
    % NOTE: When renaming the class name Untitled, the file name
    % and constructor name must be updated to use the class name.
    %
    % This template includes most, but not all, possible properties,
    % attributes, and methods that you can implement for a System object.

    % Public, nontunable properties
    properties (Nontunable)       
        tEnd = 3600                 % (s) simulation end time              
        
        Hp = 7                      % (m) height of prototype bin
        bp = 0.3214                 % inner nondimensional radius of prototype bin
        a = 0.044                    % inner nondimensional radius at top of center channel
        a0 = 0.044                   % inner nondimensional radius at outlet
        h = 0.0011                   % nondimensional height of top flow boundary        
        b = 0.3214                  % inner nondimensional radius
        H_ = 0.1                   % (m) bin height used for numerical computation
        
        thetaA = 0.5                  % ambient temperature inside
        T0 = 800                    % (°C) initial temperature of particles in bin
        Tinf = 20                   % (°C) ambient temperature
        Tref = 0                    % (°C) reference temperature 
        hInf = 10                   % (W/m2K) heat transfer coefficient to surroundings
        hcw = 10                     % (W/m2K) wall-particle boundary contact coefficient
        hcwA = 10                    % (W/m2K) wall-tank air convection coefficient        
        tauW1 = 1.7057e-04          % nondimensional time constant for wall boundary ramp function
        hp1 = 5                     % (W/m2K) top flow surface h
        hp2 = 10                    % (W/m2K) bottom overall convection coefficient for prototype bin
        hp3 = 4                     % (W/m2K) center flow channel h
        hp4 = 10                    % (W/m2K) wall overall convection coefficient for prototype bin
        hp5 = 10                    % (W/m2K) free surface h
        hp5D = 0.2                  % (W/m2K) free surface h for discharge mode
        hp5H = 0.2                  % (W/m2K) free surface h for holding mode
        hp5C = 0.1                  % (W/m2K) free surface h for chargeing mode
        hpTop = 5                   % (W/m2K) overall heat transfer coefficient for the top of the bin        
        
        kp = 0.4                    % (W/mK) particle packed thermal conductivity
        rhopPack = 2000             % (kg/m3) particle packed bulk density
        rhopLoose = 1810            % (kg/m3) particle loose bulk density
        cpp = 1025.965              % (J/kgK) average particle heat capacity  
        mup = 2.5*1.81e-5           % (kg/ms) viscosity of particles moving in air        
        
        g = 9.80665                 % (m/s2) gravitational acceleration                     
        deltaM = 0.005               % nondimensional mixing depth for charging mode        
                                                  
        clim = 1e-7         % only fourier coefficients > clim are used
        p = 4000            % total number of beta values computed
        pf = 2500           % range for beta values to be computed
        q = 4000            % total number of eta values computed
        qf = 2500           % range for eta values to be computed
        ni = 30             % number of beta values used in computation of Cnm
        mi = 40             % number of eta values used in computation of Cnm
        climFD = 5e-7       % only fourier coefficients > clim are used
        pFD = 4000          % total number of beta values computed
        pfFD = 2500         % range for beta values to be computed
        qFD = 4000          % total number of eta values computed
        qfFD = 2500         % range for eta values to be computed
        niFD = 30           % number of beta values used in computation of Cnm
        miFD = 30           % number of eta values used in computation of Cnm
        climH = 1e-6        % only fourier coefficients > clim are used
        pH = 4000           % total number of beta values computed
        pfH = 2500          % range for beta values to be computed
        qH = 2500           % total number of eta values computed
        qfH = 1000          % range for eta values to be computed
        niH = 25            % number of beta values used in computation of Cnm
        miH = 20            % number of eta values used in computation of Cnm
        climFH = 5e-7       % only fourier coefficients > clim are used
        pFH_ = 1000         % total number of beta values computed
        pfFH = 100          % range for beta values to be computed
        niFH = 25           % number of beta values used in computation of Cnm
        climW = 1e-7        % only fourier coefficients > clim are used
        climWBC = 5e-6      % only fourier coefficients > clim are used
        pW = 2000           % total number of beta values computed
        pWBC = 2000         % total number of beta values computed
        pfW = 1000          % range for beta values to be computed
        pfWBC = 1000        % range for beta values to be computed
        qW = 2000           % total number of eta values computed
        qWBC = 2000         % total number of eta values computed
        qfW = 1000          % range for eta values to be computed
        qfWBC = 1000        % range for eta values to be computed
        niW = 50            % number of beta values used in computation of Cnm
        niWBC = 50          % number of beta values used in computation of Cnm
        niWBCtot = 10       % number of beta values used in summation
        miW = 50            % number of eta values used in computation of Cnm
        objFolder = ''      % folder to save FF object               
        miWBC = 50          % number of eta values used in computation of Cnm
        miWBCtot = 10       % number of eta values used in summation
        
        ls = 40             % max storage size (time steps) for data storage file
        thetaFolder = 'thetaTest'    % folder to save theta matrices in
    end

    % properties that shouldn't be set by user
    properties (Access = protected)
        dt = 1                      % (s) time-step
        iteration = 0;
        IC                      % initial condition storage
        zIC                     % z coordinates for initial condition storage
        rIC                     % radial coordinates for initial condition storage
        rtop = []               % r-coordinates for calculations in top boundary
        rbar = []               % r-coordinates for calculations in stagnant region and top boundary
        rhat = []               % r-coordinate for calculations in center channel
        rbarH = []              % r-coordinates for 'H' and 'C' computations
        zcenter = []            % z-coordinates for calculations in center channel        
        zbar = []               % z-coordinates for calculations in stagnant region
        zbar0 = []              % initial z-coordinates for stagnant region
        zhat = []               % z-coordinates for calculations in the top boundary       
        zbarH = []              % z-coordinates for 'H' and 'C' computationses
        nzIC = 2000             % number of z-nodes for full z coordinate storage
        nrIC = 1000             % number of r-nodes for full r coordinate storage
        nrbar = 600             % number of r-nodes to use for stagnant region
        nrH = 600               % number of r-nodes to use for 'H' and 'C' computations       
        nzbar                   % number of z-nodes to use for stagnant region 
        nzbar0 = 1000           % max number of z-nodes for stagnant region
        nzH                     % number of z-nodes to use for 'H' and 'C' computations
        nzH0 = 2000             % max number of z-nodes to use for 'H' and 'C' modes
        nzc                     % number of nodes to use for center flow channel
        nzc0 = 100              % max number of z-nodes for center channel
        nzhat = 25              % number of z-nodes to use for top flow channel
        nrhat = 50              % number of r-nodes for center channel
        nzbarW = 1000           % number of nodes used in zbarW
        nrbarW = {100}          % number of nodes used in rbarW
        z = []                  % z-coordinates for whole domain
        r = []                  % r-coordinates for whole domain
        zmc = []                % z-mesh for mass flow cone
        dr = 0.01                  % radial mesh size for full domain computations
        dz = 0.01                  % z-mesh size for full domain computations
        drH = 0.001                 % nondimensional radius mesh for 'H'  
        drtop = 0.005                % nondimensional radius top boundary mesh size
        drbar = 0.00005               % nondimensional radius large mesh size
        drhat = 0.002                % nondimensional radius small mesh size
        dzH = 0.001                 % nondimensional height mesh for 'H'
        dzc = 0.005                  % nondimensional height center channel mesh size
        dzbar = 0.0001               % nondimensional height large mesh size         
        dzhat = 0.0001                % nondimensional height small mesh size
        Fo = []                 % non-dimensional time (Fourier number) vector
        FoEnd = 0               % nondimensional end time
        FoNow = 0              % current iteration non-dimensional time
        FoModePrev = 'H'          % current cycle mode (H, D, C)
        FoEmpty = 0            % non-dimensional time when tank is completely empty
        tEmpty = 0             % time when tank is completely empty
        beta = []               % stagnant radial temperature solution
        eta = []                % stagnant z temperature solution
        etaFD = []              % eigenvalues for discharge filter 1
        betaFD = []             % eigenvalues for discharge filter 2
        betaH = []              % eigenvalues for holding temperature solution
        etaH = []               % eigenvalues for holding temperature solution
        betaFH = []             % eigenvalues for holding solution filter 1
        betaC = []              % eigenvalues for charging solution
        etaC = []               % eigenvalues for charging solution
        betaFC = []             % eigenvalues for charging solution filter 1
        PsiCnm = []             % homogeneous temperature solution coefficients
        GCnm = []               % Green's Function coefficients
        FDCn = []               % Coefficients for stagnant solution filter 1
        FDCm = []               % Coefficients for stagnant solution filter 2
        CnmH = []               % Coefficients for holding solution
        FHhatCn = []            % Coefficients for holding solution filter 1
        CnmC = []               % Coefficients for charging solution
        FCCn = []               % Coefficients for charging solution filter 1      
        psiS = []               % stagnant region homogeneous temperature solution         
        thetaO = []             % centerline temperature at outlet
        thetaOB = []            % bulk temperature at outlet
        thetaS = []                  % stagnant region temperature solution
        thetaSDot = []          % stagnant region temperature time derivative
        FS0 = []                % initial condition matrix for stagnant region
        FD = []                 % discharg homogeneaous filtering function
        thetaIC = []            % initial condition contribution in Green's Formula
        thetaICDot = []         % IC contribution derivative
        RIC = []                % radial initial condition component (static)
        thetaBC1 = []           % stagnant temperature contribution from BC1
        thetaBC1Dot = []        % derivative of stagnant temperature contribution from BC1
        g1 = []                  % nonhomogeneous time-dependent BC1
        g2 = []                  % nonhomogeneous BC2
        thetaBC3 = []           % stagnant temperature contribution from BC3
        thetaBC3Dot = []        % derivative of stagnant temperature contribution from BC3
        g3 = []                  % nonhomogeneous time-dependent BC3
        g4 = []                  % nonhomogeneous BC4
        g4f = []                 % s.s. filter portion of BC4
        thetaT = []             % top boundary temperature solution
        expTop = []             % static matrix exponential for top solution
        expMC = []              % static matrix exponential for mass flow cone
        topLoss = []            % temperature loss accross top moving surface
        FT0 = []                % initial temperature for top boundary
        thetaC = []             % center boundary temperature solution
        scPe = []               % sensitivity w.r.t. Pe
        scQc = []               % sensitivity w.r.t. Qc
        sca = []                % sensitivity w.r.t. a
        FC0 = []                % initial temperature for center boundary
        thetaChat = []          % boundary condition at zhat = 1 for center channel
        thetaI = []             % offset boundary temperature for stagnant region
        thetaWI = []            % offset boundary temerature for wall transient
        rhoH = []               % initial condition for holding solution
        rhoMC = []              % initial condition for cone temperature        
        FH = []                 % filtering function for holding solution
        FCh = []                % filtering function for charging solution
        thetaH = []             % holding temperature solution sequence
        thetaCh = []            % charging temperature solution sequence
        thetaCi = []            % charging inlet particle temperature
        thetaCp = []            % prescribed temperature boundary temp for charging
        thetaMC = []            % charging mixing depth temperature
        theta = []              % full temperature solution
        thetaK = []             % full temperature solution for Kevin's Model Sim
        thetaMatchD = []        % discharge temperature with matched mesh
        thetaMatchH = []        % charge and holding temperature with matched mesh
        ubar = []               % top boundary radial velocity solution
        uzc = []                % velocity in cone for mass flow 
        wbar = []               % center boundary vertical velocity solution
        ur = []                 % radial velocity component at every cell in bin
        uz = []                 % z velocity component at every cell in bin
        Uinf = 0.02             % bulk velocity at outlet 
        AT = []                 % state matrix for top boundary temperature
        OmegaT1 = []            % diagonal elements of state matrix for top flow channel
        OmegaT2 = []            % diagonal elements of state matrix for top flow channel
        OmegaT3 = []            % diagonal elements of state matrix for top flow channel
        OmegaT4 = []            % diagonal elements of state matrix for top flow channel
        OmegaT5 = []            % diagonal elements of state matrix for top flow channel
        LambdaT1 = []           % boundary elements of state matrix for top flow channel
        LambdaT2 = []           % boundary elements of state matrix for top flow channel
        LambdaT3 = []           % boundary elements of state matrix for top flow channel
        GT = []                 % boundary state matrix for top flow channel
        AC = []                 % state matrix for center channel temperature
        AMC = []                % state matrix for mass flow cone
        dAC_Pe = []             % state partial w.r.t. Pe
        dAC_Qc = []             % state partial w.r.t. Qc
        dAC_a = []              % state partial w.r.t. a
        OmegaC1 = []            % diagonal elements of state matrix for center flow channel
        OmegaC2 = []            % diagonal elements of state matrix for center flow channel
        OmegaC = []             % diagonal elements of state matrix for center flow channel
        OmegaC4 = []            % diagonal elements of state matrix for center flow channel
        OmegaC5 = []            % diagonal elements of state matrix for center flow channel
        LambdaC1 = []           % boundary elements of state matrix for center flow channel
        LambdaC2 = []           % boundary elements of state matrix for center flow channel
        GC = []                 % boundary state matrix for center flow channel
        Rwall = {}              % cell array storing resistance of wall layers
        Rbase = []              % cell array storing resistance of base layers
        Rtop = []               % cell array storing resistance of top layers
        CapWall = {}            % cell array storing capacitance of wall layers
        CapBase = []            % cell array storing capacitance of base layers
        CapTop = []             % cell array storing capacitance of top layers
        thetaWall = zeros(4, 1)          % state variables for wall
        thetaBase = zeros(2, 1)          % state variables for base
        thetaTop = []           % state variables for top 
        baseInsulation         % insulation info for base
        wallInsulation         % insulation info for wall
        roofInsulation         % insulation info for roof of bin
        tauWall = NaN*ones(1, 2)            % non-dimensional time constants for wall layers
        tauBase = []            % non-dimensional time constants for base layers
        tauTop = []             % non-dimensional time constants for top layers
        Awall = NaN*ones(1, 1)              % state matrix for wall system
        Abase = []              % state matrix for base system
        Atop = []               % state matrix for top system
        Bwall = NaN*ones(1, 1)              % input matrix for wall system
        Bbase = []              % input matrix for base system
        Btop = []               % input matrix for top system
        Cwall = NaN*ones(2, 1)              % state output matrix for wall system
        Cbase = []              % state output matrix for base system
        Ctop = []               % state output matrix for top system
        Dwall = NaN*ones(2, 1)              % input output matrix for wall system
        Dbase = []              % input output matrix for base system
        Dtop = []               % input output matrix for top system
        wallSys = ss(ones(1, 1), ones(1, 1), ones(2, 1), ones(2, 1))            % wall state-space system
        baseSys = []            % base state-space system
        topSys = []             % top state-space system
        thetaW = []             % composite wall solution/s        
        thetaWIC = []           % IC contribution for thetaW
        thetaWBC = []           % BC contribution for thetaW
        gW1 = 0                % current inner wall boundary condition
        IBCW = []               % stored inner wall boundary conditions, f(z, t)
        rhoW = []               % composite wall initial condition/s
        zbarW = []              % z-dimension vectors for composite wall
        rbarW = []              % r-dimension vectors for composite wall        
        AWm = []                % r-BVP boundary condition matrix
        bWm = []                % r eigenfunction coefficients (IC contribution)
        bWBCm = []              % r eigenfunction coefficients (BC contribution)
        cWr = []                % lhs of reduced r-BVP boundary condition system
        etaW = []               % r-BVP eigenvalues reduced for IC contribution
        betaW = []              % z-BVP eigenvalues reduced for IC contribution
        etaWBC = []             % r-BVP eigenvalues reduced for BC contribution
        betaWBC = []            % z-BVP eigenvalues reduced for BC contribution
        Biw1 = []               % biot number for inner-most composite layer
        Biw1A  = []             % "" for tank air connection
        Biw2 = []               % biot number for outer-most composite layer        
        betaWStatic = []        % static betaW values
        etaWStatic = []         % static etaW values
        AWmStatic = []          % boundary condition matrix for etaWStatic
        bWmStatic = []          % r-BVP coefficients for etaWStatic
        GWCnm = []              % wall coefficients for IC contribution
        GWCm = []               % 1D wall coefficients for IC contribution
        GWBCCnm = []            % wall coefficients for BC contribution
        continuity = []         % cell-to-cell evaluation of 2D continuity equation
        energy = []             % cell-to-cell evaluation of 2D energy equation
        RStatic = []            % r-dimension eigenfunctions for all eta values
        RFDStatic = []          % r-dimension discharge filter eigenfunctions for all eta values
        RC3Static = []          % r-dimension eigenfunctions used for BC3
        RC4Static = []          % r-dimension eigenfunctions used for BC4
        etaStatic = []          % static eta eigenvalues
        etaFDStatic = []        % static eta eigenvalues for discharge filter
        etaHStatic = []         % static eta eigenvalues for holding solution
        RD = []                 % r-dimension eigenfunctions for current eta values
        RFD = []                % r-dimension eigenfunctions for discharge filter with current eta values
        RC3 = []                % r-dimension eigenfunctions for boundary contribution for current eta values
        RC4 = []                % r-dimension eigenfunctions for boundary contribution for current eta values
        Qc = 0                  % center channel flow rate
        Qt = 0                  % top boundary flow rate
        qTopP = 0               % (W/m2) heat loss from top of particle region
        qWall = 0               % (W/m2) heat loss from outer-most composite layer
        qLossW = 0              % (W/m2) prescribed heat flux at bin wall
        qLossB = 0              % (W/m2) prescribed heat flux at bin base
        qLossT = 0              % (W/m2) prescribed heat flux at bin top
        simMass = []            % actual mass contained in the bin at each time step
        thryMass = []           % theoretical mass contained at each time step
        zbarExp_ = []           % matched experimental data height locations
        zbarExp = []            % experimental data height locations
        rbarExp_ = []           % matched experimental data radial locations
        rbarExp = []            % experimental data radial locations
        FoExp_ = []             % matched experimental time data
        FoExp = []              % experimental time data
        thetaExp = []           % experimental temperature data
        thetaOExp = []          % centerline temperature at outlet
        thetaOBExp = []         % bulk temperature at outlet
        Bip1 = 0                % Biot number at boundary 1
        Bip2 = 0                % Biot number at boundary 2
        Bip3 = 0                % Biot number at boundary 3
        Bip4 = 0                % Biot number at boundary 4
        Bip5 = 0                % Biot number at surface of top boundary open to air
        Bip5D = 0               % Biot number for discharge mode
        Bip5H = 0               % Biot number for holding mode
        Bip5C = 0               % Biot number for chargeing mode
        rhoPack = 2000              % (kg/m3) particle packed bulk density
        rhoLoose = 1810             % (kg/m3) particle loose bulk density
        cp = 1025.965               % (J/kgK) average particle heat capacity       
        k = 0.4                     % (W/mK) particle packed thermal conductivity
        Gap = 0                 % Galilei number
        Prp = 0                 % Prandtl number
        Rep = 0                 % Reynolds number w.r.t. Uinf
        Pep = 0                 % Peclet number w.r.t. Uinf 
        Uinfp = 0               % prototype outlet velocity
        Qp = 1e-6                   % (m3/s) discharge flow rate
        QChp = 1e-6                 % (m3/s) charge flow rate  
        Q = 0                   % (m3/s) prototype volumetric flow rate
        QCh = 0                 % (m3/s) prototype charging flow rate
        ztop = []               % zbar element of top surface
        z0C = 0                 % starting height for full charge
        HNow = 0                % (m) current height of top surface
        HNowC = 0               % (m) current height of top surface
        dH = 0                  % (m) change in height per time step
        dHCh = 0                % (m) change in height per time step for charging mode
        modZH = 0
        modZS = 0
        modZC = 0
        alphapPacked = 0        % (m2/s) particle thermal diffusivity for prototype
        alphapLoose = 0         % (m2/s) particle thermal diffusivity for prototype 
        alphaPacked = 0         % (m2/s) particle thermal diffusivity 
        alphaLoose = 0          % (m2/s) particle thermal diffusivity
        nup = 0                 % (m2/s) dynamic viscosity of particles moving in air 
        nu = 0                  % (m2/s) dynamic viscosity of particles moving in air 
        CapS = 0               % (J/m2K) relative thermal capacitance of stagnant region
        epsilon2 = 0           % ratio of base thermal capacitance to stagnant region thermal capacitance
        epsilon4 = 0           % ratio of wall thermal capacitance to stagnant region thermal capacitance 
        h1 = 10                     % (W/m2K) heat transfer coefficient at boundary 1
        h2 = 0.07                   % (W/m2K) heat transfer coefficient at boundary 2
        h3 = 5                      % (W/m2K) heat transfer coefficient at boundary 3
        h4 = 0.07                   % (W/m2K) heat transfer coefficient at boundary 4
        h5 = 0.2                     % (W/m2K) heat transfer coefficient for free surface
        h5D = 0.2                    % (W/m2K) heat transfer coefficient at free surface for discharge mode
        h5H = 0.2                    % (W/m2K) heat transfer coefficient at free surface for holding mode
        h5C = 0.1                     % (W/m2K) heat transfer coefficient at free surface for chargeing mode
        Bi1 = 0                % Biot number at boundary 1
        Bi2 = 0                % "" boundary 2
        Bi3 = 0                % "" boundary 3
        Bi4 = 0                % "" boundary 4
        Bi5 = 0                % "" surface of top boundary open to air
        Bi5D = 0               % "" for discharge mode
        Bi5H = 0               % "" for holding mode
        Bi5C = 0               % "" for chargeing mode
        Ga = 0                 % Galilei number
        Pr = 0                 % Prandtl number
        Re = 0                 % Reynolds number w.r.t. Uinf
        Pe = 0                 % Peclet number w.r.t. Uinf 
        cGet = []               % index array for obtaining eigenvalues
        cGetFD = []             % index array for obtaining eigenvalues
        cGetH = []              % index array for obtaining eigenvalues
        cGetFH = []             % index array for obtaining eigenvalues
        cGetW = []              % index array for obtaining eigenvalues
        cGetWBC = []            % index array for obtaining eigenvalues
        df = 0                 % nondimensional time-step
        vtbl = []               % table containing all static variables
        cfig = []               % figure showing Fourier coefficients
        betafig = []            % figure showing beta values
        etafig = []             % figure showing eta values
        psifig = []             % figure showing temperature for psi
        thetafig = []           % figure showing temperature for theta (discharge)
        thetaHfig = []          % figure showing temperature for theta (holding)
        thetaCfig = []          % figure showing temperature for theta (charging)
        ubarfig = []            % figure showing velocity in top boundary
        ubarContFig = []        % figure showing continuity solution in top bound.
        wbarfig = []            % figure showing velocity in center boundary
        computeBUfig = []       % figure showing approximated bulk velocities
        thetaOfig = []          % figure showing bulk outlet temperature
        betaHfig = []           % figure showing holding beta eigenvalues
        etaHfig = []            % figure showing holding eta eigenvalues
        etaFHfig = []           % figure showing holding filter eta eigenvalues
        cfigFHhat = []          % figure showing holding filter coefficients
        cfigH = []              % figure showing holding Fourier coefficients
        betaCfig = []           % figure showing charging beta eigenvalues
        etaCfig = []            % figure showing charging eta eigenvalues
        etaFCfig = []           % figure showing charging filter eta eigenvalues
        cfigFC = []             % figure showing charging filter coefficients
        cfigC = []              % figure showing charging Fourier coefficients
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
          icon = sprintf('Particle\nTES\nBin'); 
       end
%         function icon = getIconImpl(~)
%             % Define icon for System block
%             icon = matlab.system.display.Icon('png-transparent-simulink-matlab-mathworks-computer-software-logo-coder-miscellaneous-angle-rectangle.png');
%         end
        function [in1name, in2name, in3name] = getInputNamesImpl(~)
          in1name = 'Tin';
          in2name = 'mdot';
          in3name = 't';
        end
        function [out1name, out2name, out3name, out4name] = getOutputNamesImpl(~)
          out1name = 'Tout';
          out2name = 'Tbulk';
          out3name = 'Estored';
          out4name = 'ztop';
        end             
        function groups = getPropertyGroupsImpl
          group1 = matlab.system.display.SectionGroup( ...
              'Title', 'Mesh Parameters', ...
              'PropertyList', {});          
          group2 = matlab.system.display.SectionGroup( ...
              'Title', 'Bin Geometry Parameters', ...
              'PropertyList', {'Hp', 'bp', 'a', 'h'});
          group3 = matlab.system.display.SectionGroup( ...
              'Title', 'Heat Transfer Parameters', ...
              'PropertyList', {'T0', 'Tinf', 'Tref', 'hInf', 'hcw', 'hcwA', ...
              'tauW1', 'hp1', 'hp2', 'hp3', 'hp4', 'hp5', 'hp5D', 'hp5H', ...
              'hp5C', 'hpTop'});
          group4 = matlab.system.display.SectionGroup( ...
              'Title', 'Particle Properties', ...
              'PropertyList', {'kp', 'rhopPack', 'rhopLoose', 'cpp', 'mup'});                             
          groups = [group2, group3, group4];    
       end
       
    end

    methods(Access = protected)
        %% Common functions
        function setupImpl(obj)
            % Perform one-time calculations, such as computing constants
            obj.dt = 1;
            obj.FoEnd = t2Fo(obj, obj.tEnd);
            obj.df = t2Fo(obj, obj.dt);
%             initializeWallSys(obj);
            obj.ztop = 0.05;
            obj.rtop = obj.a:obj.drtop:obj.b;
            [obj.rbar, obj.drbar] = ... 
                        nodeGen(obj, [obj.a0, obj.b], obj.nrbar);
            [obj.rbarH, obj.drH] = ... 
                        nodeGen(obj, [0, obj.b], obj.nrH);
            obj.rhat = 1e-6:obj.drhat:obj.a0;
            obj.r = [obj.rhat, obj.rbar(2:end)];
            obj.nzbar = ceil(obj.nzbar0*obj.ztop);
            obj.zbar0 = nodeGen(obj, [0, obj.ztop], obj.nzbar);
            [obj.zbar, obj.dzbar] = ... 
                        nodeGen(obj, [0, obj.ztop], obj.nzbar);
            obj.nzH = ceil(obj.nzH0*obj.ztop);
            [obj.zbarH, obj.dzH] = ... 
                        nodeGen(obj, [0, obj.ztop], obj.nzH);
            obj.zcenter = 0:obj.dzc:obj.ztop;
            obj.zhat = 0:obj.dzhat:obj.h;
            obj.z = [obj.zbar0, 1+obj.zhat(2:end)];
            obj.zIC = NaN*ones(1, obj.nzIC);
            obj.rIC = NaN*ones(1, obj.nrIC);
            obj.IC = NaN*ones(obj.nzIC, obj.nrIC);
            obj.thetaS = ones(length(obj.zbar), length(obj.rbar));
            obj.HNow = obj.zbar(end)*obj.H_;    
            obj.HNowC = obj.zbarH(end)*obj.H_;
            obj.alphapPacked = obj.kp/(obj.rhopPack*obj.cpp); 
            obj.alphapLoose = obj.kp/(obj.rhopLoose*obj.cpp);
            obj.alphaPacked = obj.k/(obj.rhoPack*obj.cp); 
            obj.alphaLoose = obj.k/(obj.rhoLoose*obj.cp);
            obj.CapS = obj.rhoPack*obj.cp*obj.H_;
%             obj.epsilon2 = obj.CapBase/obj.CapS;
%             obj.epsilon4 = obj.CapWall/obj.CapS;
            obj.nup = obj.mup/obj.rhoLoose;
            obj.nu = obj.nup;
            obj.Bi1 = obj.h1*obj.H_/obj.k;          
            obj.Bi2 = obj.h2*obj.H_/obj.k;               
            obj.Bi3 = obj.h3*obj.H_/obj.k;               
            obj.Bi4 = obj.h4*obj.H_/obj.k;                    
            obj.Bi5D = obj.h5D*obj.H_/obj.k;  
            obj.Bi5H = obj.h5H*obj.H_/obj.k; 
            obj.Bi5C = obj.h5C*obj.H_/obj.k; 
            obj.Ga = obj.g*obj.H_^3/obj.nu^2;                
            obj.Pr = obj.nu/obj.alphaPacked;                        
%             obj.Uinf = obj.Q/(pi*(obj.a0*obj.H_)^2);
%             obj.Q = obj.Uinf*pi*(obj.a0*obj.H_)^2;
%             obj.Re = obj.H_*obj.Uinf/obj.nu;           
%             obj.Pe = obj.Re*obj.Pr;
%             obj.Qc = obj.Q/(obj.H_^2*obj.Uinf);
            obj.Bip1 = obj.hp1*obj.Hp/obj.kp;          
            obj.Bip2 = obj.hp2*obj.Hp/obj.kp;               
            obj.Bip3 = obj.hp3*obj.Hp/obj.kp;               
            obj.Bip4 = obj.hp4*obj.Hp/obj.kp;                    
            obj.Bip5D = obj.hp5D*obj.Hp/obj.k;  
            obj.Bip5H = obj.hp5H*obj.Hp/obj.k; 
            obj.Bip5C = obj.hp5C*obj.Hp/obj.k;   
            obj.Gap = obj.g*obj.Hp^3/obj.nup^2;                
            obj.Prp = obj.nup/obj.alphapPacked;                        
            obj.Uinfp = obj.Qp/(pi*(obj.a0*obj.Hp)^2);
            obj.Rep = obj.Hp*obj.Uinfp/obj.nup;           
            obj.Pep = obj.Rep*obj.Prp;
            obj.df = obj.t2Fo(obj.dt);
            obj.dH = obj.Q*obj.dt/(pi*(obj.b*obj.H_)^2);
            obj.dHCh = obj.QCh*obj.dt/(pi*(obj.b*obj.H_)^2);
            obj.modZH = min(ceil(obj.dzH./(obj.dH/obj.H_)));
            obj.modZS = min(ceil(obj.dzbar./(obj.dH/obj.H_)));
            obj.modZC = min(ceil(obj.dzc./(obj.dH/obj.H_)));
            obj.tEmpty = obj.H_*pi*(obj.b*obj.H_)^2/obj.Q;
            obj.FoEmpty = obj.t2Fo(obj.tEmpty);
            obj.Fo = 0:obj.df:obj.FoEnd;
            obj.FoNow = 0;
            % initialization functions
            setInsulation(obj);
            matchSimilarityParams(obj);
        end
        function resetImpl(obj)
            % Initialize / reset discrete-state properties
            % recomputes static variables. needs to be run if any of the
            % system properties are changed externally.
            obj.rtop = obj.a:obj.drtop:obj.b;
            [obj.rbar, obj.drbar] = ... 
                        nodeGen(obj, [obj.a0, obj.b], obj.nrbar);
            [obj.rbarH, obj.drH] = ... 
                        nodeGen(obj, [0, obj.b], obj.nrH);
            obj.rhat = 1e-6:obj.drhat:obj.a0; 
            obj.r = [obj.rhat, obj.rbar(2:end)];
            obj.zbar0 = nodeGen(obj, [0, 1], obj.nzbar);
            [obj.zbar, obj.dzbar] = ... 
                        nodeGen(obj, [0, obj.ztop], obj.nzbar);
            obj.zcenter = 0:obj.dzc:obj.ztop;
            obj.zhat = 0:obj.dzhat:obj.h; 
            [obj.zbarH, obj.dzH] = ... 
                        nodeGen(obj, [0, obj.ztop], obj.nzH);
            obj.z = [obj.zbar0, 1+obj.zhat(2:end)];
            obj.HNow = obj.zbar(end)*obj.H_;  
            obj.HNowC = obj.zbarH(end)*obj.H_;
            obj.alphapPacked = obj.kp/(obj.rhopPack*obj.cpp); 
            obj.alphapLoose = obj.kp/(obj.rhopLoose*obj.cpp); 
            obj.alphaPacked = obj.k/(obj.rhoPack*obj.cp); 
            obj.alphaLoose = obj.k/(obj.rhoLoose*obj.cp);  
            obj.nup = obj.mup/obj.rhoLoose;
            obj.CapS = obj.rhoPack*obj.cp*obj.H_;
%             obj.epsilon2 = obj.CapBase/obj.CapS;
%             obj.epsilon4 = obj.CapWall/obj.CapS;
            obj.Bi1 = obj.h1*obj.H_/obj.k;          
            obj.Bi2 = obj.h2*obj.H_/obj.k;               
            obj.Bi3 = obj.h3*obj.H_/obj.k;               
            obj.Bi4 = obj.h4*obj.H_/obj.k;                    
            obj.Bi5D = obj.h5D*obj.H_/obj.k;  
            obj.Bi5H = obj.h5H*obj.H_/obj.k; 
            obj.Bi5C = obj.h5C*obj.H_/obj.k; 
            obj.Ga = obj.g*obj.H_^3/obj.nu^2;                
            obj.Pr = obj.nu/obj.alphaPacked;   
%             obj.Uinf = obj.Q/(pi*(obj.a0*obj.H_)^2);
%             obj.Q = obj.Uinf*pi*(obj.a0*obj.H_)^2;
            obj.Qc = obj.Q/(obj.H_^2*obj.Uinf);
            obj.Re = obj.H_*obj.Uinf/obj.nu;                 
            obj.Pe = obj.Re*obj.Pr;
            obj.Bip1 = obj.hp1*obj.Hp/obj.kp;          
            obj.Bip2 = obj.hp2*obj.Hp/obj.kp;               
            obj.Bip3 = obj.hp3*obj.Hp/obj.kp;               
            obj.Bip4 = obj.hp4*obj.Hp/obj.kp;                    
            obj.Bip5D = obj.hp5D*obj.Hp/obj.k;  
            obj.Bip5H = obj.hp5H*obj.Hp/obj.k; 
            obj.Bip5C = obj.hp5C*obj.Hp/obj.k;  
            obj.Gap = obj.g*obj.Hp^3/obj.nup^2;                
            obj.Prp = obj.nup/obj.alphapPacked;                        
            obj.Uinfp = obj.Qp/(pi*(obj.a0*obj.Hp)^2);
            obj.Rep = obj.Hp*obj.Uinfp/obj.nup;           
            obj.Pep = obj.Rep*obj.Prp;
            obj.dH = obj.Q*obj.dt/(pi*(obj.b*obj.H_)^2);
            obj.dHCh = obj.QCh*obj.dt/(pi*(obj.b*obj.H_)^2);
            obj.modZH = min(ceil(obj.dzH./(obj.dH/obj.H_)));
            obj.modZS = min(ceil(obj.dzbar./(obj.dH/obj.H_)));
            obj.modZC = min(ceil(obj.dzc./(obj.dH/obj.H_)));
            obj.tEmpty = obj.ztop*obj.H_*pi*(obj.b*obj.H_)^2/obj.Q;
            obj.FoEmpty = obj.t2Fo(obj.tEmpty);
            obj.Fo = 0:obj.df:obj.FoEnd;
        end
        function [y1, y2, y3, y4] = stepImpl(obj, Tin, mdot, t)
            % Implement algorithm. Calculate y as a function of input u and
            % discrete states.
            % initialize temperature matrices
            obj.df = obj.t2Fo(t, 1) - obj.FoNow;
            obj.dt = obj.Fo2t(obj.df, 0);             
            if obj.FoNow == 0
                obj.iteration = 1;
                if mdot > 0
                    obj.ztop = 0.1;
                else
                    obj.ztop = 0.9;
                end
                [r_, obj.drH] = nodeGen(obj, [0, obj.b], obj.nrH);
                obj.nzH = ceil(obj.nzH0*obj.ztop);
                [z_, obj.dzH] = nodeGen(obj, [0, obj.ztop], obj.nzH);
                nch = length(z_); mch = length(r_);
                IC_ = ones(nch, mch);
                initializeBaseSys(obj);
                initializeWallSys(obj);
                computeBaseSys(obj, 1, obj.df);
                computeWallSys(obj, 1, obj.df);
            else
                [i, j] = find(~isnan(obj.IC)); n = max(i); m = max(j);                
                IC_ = obj.IC(1:n, 1:m);
                z_ = obj.zIC(1:n);
                r_ = obj.rIC(1:m);
            end          
            % set cycle mode
            if mdot == 0, FoMode = 'H'; end
            if mdot > 0, FoMode = 'C'; end
            if mdot < 0, FoMode = 'D'; end
            % set limits on filling and emptying
            if obj.ztop < 0.05
                if FoMode == 'D'
                    FoMode = 'H';
                end                    
            end
            if obj.ztop > 0.95
                if FoMode == 'C'
                    FoMode = 'H';
                end                    
            end
            % simulate cycle time step
            switch FoMode
                case 'H'
                    % run holding simulation
                    if obj.FoModePrev ~= 'H'
                        matchHoldIC(obj, IC_, z_, r_);
                    end
                    obj.Bi5 = obj.Bi5H;
                    % compute hold temperature values for current time step                                                                               
                    theta_ = computeThetaH(obj, obj.df, IC_, z_, r_);                
                    % update initial condition
                    IC_ = theta_;                                               
                case 'C'
                    % run charging simulation
                    if obj.FoModePrev ~= 'C'
                        [IC_, z_, r_] = matchChargeIC(obj, IC_, z_, r_);
                    end
                    % set charge flow rate parameters
                    obj.QChp = abs(mdot)/obj.rhopPack;
                    Uinfp_ = obj.QChp/(pi*(obj.a0*obj.Hp)^2);
                    Uinf_ = obj.Hp/obj.H_*Uinfp_;
                    obj.QCh = Uinf_*(pi*(obj.a0*obj.H_)^2);
                    obj.thetaCi = (Tin - obj.Tinf)/(obj.T0 - obj.Tinf);
%                     resetImpl(obj);
                    obj.Bi5 = obj.Bi5C;
                    % compute charge temperature values for current time step
                    theta_ = computeThetaH(obj, obj.df, IC_, z_, r_);
                    % update spatial domain
                    obj.HNow = obj.HNow + obj.dHCh;
                    obj.ztop = obj.HNow/obj.H_;
                    obj.nzH = ceil(obj.nzH0*obj.ztop);
                    [z_, obj.dzH] = nodeGen(obj, [0, obj.ztop], obj.nzH);
                    % update temperature domain
                    move = abs(size(theta_, 1) - length(z_));
                    if move > 0             
                        % update temperature 
%                         for n = 1:move
                            theta_ = [theta_; ...
                                   obj.thetaCi*ones(move, length(theta_(1, :)))];
                            % compute temperature distribution for mixing depth
%                             [~, im] = min(abs(z_ - (obj.ztop - obj.deltaM)));
%                             zm = z_(im:end);
%                             if length(zm) < 3
%                                 im = im - 3;
%                                 zm = z_(im:end);
%                             end
%                             fz = simpsonIntegrator(obj, zm);
%                             thetaMC_ = 1/(zm(end)-zm(1))*fz*theta_(im:end, :);
%                             theta_(im:end, :) = repmat(thetaMC_, length(zm), 1);
%                         end
                    end 
                    % update initial condition
                    IC_ = theta_;
                
                case 'D'
                    % run discharging simulation
                    if obj.FoModePrev ~= 'D' || obj.FoNow == 0
                        [thetaS_, zbar_, thetaC_, zcenter_, thetaT_] ...
                            = matchDischargeIC(obj, IC_, z_, r_);
                    else
                        [thetaS_, zbar_, thetaC_, zcenter_, thetaT_] ...
                             = disectDischargeDomains(obj, IC_, z_, r_);                       
                    end
                    % set discharge flow rate parameters
                    obj.Qp = abs(mdot)/obj.rhopPack;
                    obj.Uinfp = obj.Qp/(pi*(obj.a0*obj.Hp)^2);
                    matchSimilarityParams(obj);
                    obj.Bi5 = obj.Bi5D;
                    % compute temperature in top boundary
                    thetaT_ = iterateThetaT(obj, thetaS_, obj.rbar, thetaT_, obj.zhat, obj.rtop, zcenter_);
                    % compute temperature in center boundary
                    thetaC_ = iterateThetaC(obj, thetaS_, zbar_, thetaC_, zcenter_, obj.rhat);
                    % compute the stagnant boundary offset temperature
                    thetaI_ = computeThetaI(obj, thetaT_, thetaC_, obj.rtop);
                    % compute temperature in stagnant region
                    thetaS_ = iterateThetaS(obj, thetaC_, zcenter_, thetaT_, obj.rtop, thetaI_, thetaS_, zbar_, obj.rbar);                                                           
                    % update spatial domain
                    obj.HNow = obj.HNow - obj.dH;
                    obj.ztop = obj.HNow/obj.H_;
                    obj.nzbar = ceil(obj.nzbar0*obj.ztop);
                    [zbar_, obj.dzbar] = nodeGen(obj, [0, obj.ztop], obj.nzbar);
                    zcenter_ = 0:obj.dzc:obj.ztop;
                    % update temperature array
                    move = abs(size(thetaS_, 1) - length(zbar_));
                    movec = abs(size(thetaC_, 1) - length(zcenter_));
                    if move > 0             
                        % update temperature 
                        for n = 1:move
                            for i = 1:length(obj.zhat)-1
                                thetaT_(i, :) = thetaT_(i+1, :);                    
                            end
                            thetaT_(end, :) = obj.S2T(thetaS_, obj.rbar, obj.rtop);
                            thetaS_ = thetaS_(1:end-1, :); 
                        end
                    end
                    if movec > 0
                        for n = 1:movec
                            thetaC_ = thetaC_(1:end-1, :);                   
                        end
                    end                   
                    % store temperature values in solution matrix
                    theta_ = patchTheta(obj, thetaS_, thetaT_, thetaC_, zbar_, obj.rbar, obj.zhat, zcenter_, obj.rtop, obj.rhat);
                    % update initial condition
                    IC_ = theta_; 
                    z_ = [zbar_, zbar_(end) + obj.zhat(2:end)];
                    r_ = [obj.rhat(1:end-1), obj.rbar];                  
            end          
                   
            % update wall and base systems
            ub = mean(theta_(1, :));
            uw = mean(theta_(:, end));
            computeBaseSys(obj, ub, obj.df);
            computeWallSys(obj, uw, obj.df);
            
            % calculate bulk outlet temperature
            [~, ri] = min(abs(obj.a0 - r_));
            fr = simpsonIntegrator(obj, r_(1:ri));
            Tout = obj.theta2T(2/r_(ri)^2*(IC_(1, 1:ri).*r_(1:ri))*fr');
            
            % calculate bulk volumetric temperature
            fz = simpsonIntegrator(obj, z_);
            fr = simpsonIntegrator(obj, r_);
            Ir = 2*(r_.*IC_)*fr';
            Tbulk = obj.theta2T(Ir'*fz'./(max(r_)^2*max(z_)));
            
            % calculate total stored energy in bin
            Estored = obj.rhoPack*pi*max(r_)^2*max(z_)*obj.cp*Tbulk;
            
            % set outputs
            y1 = Tout;
            y2 = Tbulk;
            y3 = Estored;
            y4 = obj.ztop;
            
            % update initial condition
            obj.FoModePrev = FoMode;
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
        function flag = isInactivePropertyImpl(obj,prop)
            % Return false if property is visible based on object 
            % configuration, for the command line and System block dialog
            flag = false;
        end        
        %% model functions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % TES setup functions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function setInsulation(obj)
            % defines wall and base insulation parameters (hard coded for
            % now
            % set insulation specifications for bin wall and base
            obj.wallInsulation = cell(4, 5);
            obj.wallInsulation{1, 1} = 'tufcrete 47';
            obj.wallInsulation{1, 2} = [2.25, 2.35]; % m
            obj.wallInsulation{1, 3} = 1.53;         % W/mK
            obj.wallInsulation{1, 4} = 2210;         % kg/m3
            obj.wallInsulation{1, 5} = 1175;         % J/kgK
            obj.wallInsulation{2, 1} = 'skamolex';
            obj.wallInsulation{2, 2} = [2.35, 2.55]; % m
            obj.wallInsulation{2, 3} = 0.09;         % W/mK
            obj.wallInsulation{2, 4} = 245;          % kg/m3
            obj.wallInsulation{2, 5} = 840;          % J/kgK
            obj.wallInsulation{3, 1} = 'elmtherm';
            obj.wallInsulation{3, 2} = [2.55, 2.65];  % m
            obj.wallInsulation{3, 3} = 0.025;        % W/mK
            obj.wallInsulation{3, 4} = 270;          % kg/m3
            obj.wallInsulation{3, 5} = 1005;         % J/kgK
            obj.wallInsulation{4, 1} = 'ss304';
            obj.wallInsulation{4, 2} = [2.65, 2.65635]; % m
            obj.wallInsulation{4, 3} = 30;             % W/mK
            obj.wallInsulation{4, 4} = 7700;           % kg/m3
            obj.wallInsulation{4, 5} = 500;            % J/kgK
            % set insulation specifications for bin base
            % new insulation configuration
            obj.baseInsulation = cell(2, 5);
            obj.baseInsulation{1, 1} = 'particles';
            obj.baseInsulation{1, 2} = [0, 0.1];
            obj.baseInsulation{1, 3} = 0.4;          % W/mK
            obj.baseInsulation{1, 4} = 2000;          % kg/m3
            obj.baseInsulation{1, 5} = 1025.965;        % J/kgK
            obj.baseInsulation{2, 1} = 'fondag';
            obj.baseInsulation{2, 2} = [0.1, 0.1+0.1905];
            obj.baseInsulation{2, 3} = 1.75;          % W/mK
            obj.baseInsulation{2, 4} = 2210;          % kg/m3
            obj.baseInsulation{2, 5} = 1046.7;        % J/kgK
        end
        function matchSimilarityParams(obj)
            % matches parameters for scaled model simulation
            % Reynolds number for momentum scaling (finding optimal viscosity)
            obj.Rep = obj.Hp*obj.Uinfp/obj.nup;
            obj.nu = obj.nup;
            obj.Uinf = obj.nu*obj.Rep/obj.H_;

            % Prandlt number for energy scaling (finding optimal thermal diffusivity)
            obj.alphaPacked = obj.nu/obj.Prp;

            % thermal conductivity is relatively unaffected by particle selection at
            % high temperatures, so keep equivalent for model and prototype
            obj.k = obj.kp;
            obj.cp = obj.cpp;

            % reset model flow rates
            obj.rhoPack = obj.k/(obj.alphaPacked*obj.cp);
            obj.Q = obj.Uinf*pi*(obj.a0*obj.H_)^2;
%             obj.QCh = Qcharge/Qdischarge*obj.Q;

            % Biot numbers for boundary condition scaling
            obj.h1 = obj.Bip1*obj.k/obj.H_;
            obj.h2 = obj.Bip2*obj.k/obj.H_;
            obj.h3 = obj.Bip3*obj.k/obj.H_;
            obj.h4 = obj.Bip4*obj.k/obj.H_;
            obj.h5D = obj.Bip5D*obj.k/obj.H_;
            obj.h5H = obj.Bip5H*obj.k/obj.H_;
            obj.h5C = obj.Bip5C*obj.k/obj.H_;
            resetImpl(obj);
        end                
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % functions used to transfer between storage modes
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function matchWallBoundary(obj, z_, theta_)
            % generates the inner wall prescribed temperature BC with the
            % current air and particle temperatures
            [~, h_] = min(abs(z_(end) - obj.zbarW));
            gWs = interp1(z_, theta_, obj.zbarW(1:h_), 'makima');
%             gWs = mean(gWs);
            obj.gW1 = [gWs'.*ones(h_, 1); ...
                      obj.thetaA*ones(length(obj.zbarW) - h_, 1)];
            % add ramp function to smooth response
%             obj.gW1 = obj.gW1*(1 - exp(-obj.FoNow/obj.tauW1));
        end
        function [theta_, z_, r_] = matchHoldIC(obj, IC, zIC, rIC)
            [z_, obj.dzH] = nodeGen(obj, [0, obj.ztop], obj.nzH);
            [r_, obj.drH] =  nodeGen(obj, [0, obj.b], obj.nrH);
            [R, Z] = meshgrid(rIC, zIC);
            [Rq, Zq] = meshgrid(r_, z_);
            theta_ = interp2(R, Z, IC, Rq, Zq, 'spline');           
        end
        function [theta_, z_, r_] = matchChargeIC(obj, IC, zIC, rIC)
            [z_, obj.dzH] = nodeGen(obj, [0, obj.ztop], obj.nzH);
            [r_, obj.drH] = nodeGen(obj, [0, obj.b], obj.nrH);
            [R, Z] = meshgrid(rIC, zIC);
            [Rq, Zq] = meshgrid(r_, z_);
            theta_ = interp2(R, Z, IC, Rq, Zq, 'spline');
            theta_ = ones(size(theta_));
        end
        function [thetaS_, zbar_, thetaC_, zcenter_, thetaT_] ...
                                = matchDischargeIC(obj, IC, zIC, rIC)
            obj.ztop = obj.ztop - obj.h;
            obj.HNow = obj.ztop*obj.H_;
            obj.nzbar = ceil(obj.nzbar0*obj.ztop);
            [zbar_, obj.dzbar] = ...
                            nodeGen(obj, [0, obj.ztop], obj.nzbar);
            zcenter_ = 0:obj.dzc:obj.ztop;
            [~, i] = min(abs(zbar_(end) - zIC));
            [~, j] = min(abs(obj.rbar(1) - rIC));
            % match dimensions for stagnant region            
            [R, Z] = meshgrid(rIC(j:end), zIC(1:i));
            [Rq, Zq] = meshgrid(obj.rbar, zbar_);
            thetaS_ = interp2(R, Z, IC(1:i, j:end), Rq, Zq, 'makima');
            % match dimensions for center channel           
            [R, Z] = meshgrid(rIC(1:j), zIC(1:i));
            [Rq, Zq] = meshgrid(obj.rhat, zcenter_);
            thetaC_ = interp2(R, Z, IC(1:i, 1:j), Rq, Zq, 'makima');
            % match dimensions for top flow surface
            [R, Z] = meshgrid(rIC(j:end), zIC(i:end) - zIC(i));
            [Rq, Zq] = meshgrid(obj.rtop, obj.zhat);
            thetaT_ = interp2(R, Z, IC(i:end, j:end), Rq, Zq, 'makima');
        end
        function [thetaS_, zbar_, thetaC_, zcenter_, thetaT_] ...
                             = disectDischargeDomains(obj, IC, zIC, rIC)
            % extracts the temperature profiles and spatial domains for the
            % stagnant region, top flow surface, and center flow channel
            obj.nzbar = ceil(obj.nzbar0*obj.ztop);
            [zbar_, obj.dzbar] = ...
                            nodeGen(obj, [0, obj.ztop], obj.nzbar);
            zcenter_ = 0:obj.dzc:obj.ztop;
            [~, i] = min(abs(zbar_(end) - zIC));
            [~, j] = min(abs(obj.rbar(1) - rIC));
            % stagnant region
            thetaS_ = IC(1:i, j:end);
            % center flow channel
            [R, Z] = meshgrid(obj.rhat, zIC(1:i));
            [Rq, Zq] = meshgrid(obj.rhat, zcenter_);
            thetaC_ = interp2(R, Z, IC(1:i, 1:j), Rq, Zq, 'makima');
            % top flow channel
            [R, Z] = meshgrid(rIC(j:end), zIC(i:end) - zIC(i));
            [Rq, Zq] = meshgrid(obj.rtop, obj.zhat);
            thetaT_ = interp2(R, Z, IC(i:end, j:end), Rq, Zq, 'makima');                                   
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % conduction temperature solution for holding and charging mode
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function theta_ = computeThetaH(obj, Fo_, IC, z_, r_)
            % compute filter function solution
            FH_ = computeFH(obj, z_, r_); % needs to be before eigenvalue comp.
            % compute eigenvalues and coefficients if not yet populated
            beta_ = computeBetaH(obj);
            eta_ = computeEtaH(obj);
            [Cnm_, beta_, eta_] = computeCnmH(obj, IC, FH_, z_, r_, beta_, eta_);
            n = length(z_); m = length(r_);            
            % compute Fourier series for filtered solution
            theta_ = zeros(n, m);
            for i = 1:length(Cnm_)
                Z = XHn(obj, z_, beta_(i));
                R = XHm(obj, r_, eta_(i));
                F = Xt(obj, Fo_, beta_(i), eta_(i));
                theta_ = theta_ + Cnm_(i)*F*Z'*R;               
            end             
            % shift back with filtering function to get original solution
            theta_ = theta_ - FH_;
        end
        function FH_ = computeFH(obj, z_, r_)
            % computes the filtering function for the bin holding solution           
            if isempty(obj.g2), initializeBaseSys(obj); end
            if isempty(obj.g4), initializeWallSys(obj); end
            beta_ = computeBetaFH(obj);
            [Cn_, beta_] = computeFHhatCn(obj, beta_); 
            % compute FHhat
            n = length(z_); m = length(r_);
            % first s.s. BVP filter
            FH1 = zeros(n, m);
            for i = 1:length(Cn_)
                Z = XHn(obj, z_, beta_(i));
                R = XFHm(obj, r_, beta_(i));
                FH1 = FH1 + Cn_(i)*Z'*R;               
            end
            % second s.s. BVP filter
            FH2 = repmat(obj.g2*(z_' - obj.ztop - 1/obj.Bi5), 1, m);
            % full filter superposition
            FH_ = FH1 + FH2 - obj.thetaA;
        end
        function xi = XHn(~, zbar_, beta_)
            % eigenfunction for z-dimension BVP
            xi = cos(beta_*zbar_);
        end
        function ni = NHn(obj, beta_)
            % norm for z-dimension BVP
            ni = 0.5*(obj.ztop*(beta_.^2 + obj.Bi5^2) + obj.Bi5)./ ...
                      (beta_.^2 + obj.Bi5^2);
        end
        function zi = ZHn(obj, beta_)
            % transendental function for beta values
            zi = beta_.*sin(beta_*obj.ztop) - obj.Bi5*cos(beta_*obj.ztop);
        end 
        function xi = XHm(~, rbar_, eta_)
            % eigenfunction for r-dimension BVP
            xi = besselj(0, eta_*rbar_);
        end
        function ni = NHm(obj, eta_)
            % norm for r-dimension BVP
            if eta_ == 0
                ni = 0.5*obj.b^2;
            else
                ni = 0.5*obj.b^2*besselj(0, eta_*obj.b).^2;
            end
        end
        function zi = ZHm(obj, eta_)
            % transendental function for eta values
            zi = besselj(1, eta_*obj.b);
        end
        function xi = XFHn(obj, zbar_, eta_)
            % z-BVP solution for holding filter 2
            xi = eta_*(sinh(eta_*obj.ztop)*sinh(eta_*zbar_) - ...
                       cosh(eta_*obj.ztop)*cosh(eta_*zbar_)) + ...
                 obj.Bi5*(sinh(eta_*obj.ztop)*sinh(eta_*zbar_) - ...
                         (eta_*obj.ztop)*cosh(eta_*zbar_));
        end
        function xi = XFHm(~, rbar_, beta_)
            % r-BVP solution for holding filter 1
            xi = besseli(0, beta_*rbar_);
        end
        function c = FourierCoefficientThetaHhat(obj, IC, FH, z_, r_, beta_, eta_)
            % Fourier coefficients
            C_num = computeThetaHhatIntegral(obj, IC, FH, z_, r_, beta_, eta_);
            c = C_num/(NHm(obj, eta_)*NHn(obj, beta_));
        end
        function c = FourierCoefficientFH1(obj, beta_)
            % Fourier coefficients
            C_num = mean(mean(obj.g4))*sin(beta_*obj.ztop)./ ...
                                      (beta_.^2.*besseli(1, beta_*obj.b));
            c = C_num/NHn(obj, beta_);  
        end
        function x = computeThetaHhatIntegral(obj, IC, FH, z_, r_, beta_, eta_)
            % computes integral from numerator of Fourier coefficient for
            % the filtered thetaH solution
            fr = simpsonIntegrator(obj, r_);
            fz = simpsonIntegrator(obj, z_);
            Ir = (XHm(obj, r_, eta_).*r_.*(IC + FH))*fr';
            x = (XHn(obj, z_, beta_).*Ir')*fz';
        end                
        % computation of beta values
        function beta_ = computeBetaH(obj)
            interval = linspace(0, obj.pfH, obj.pH); % interval/spacing 
                                                     % of root calculation
            rn = NaN*ones(obj.pH, 1);                % roots initialization
            options = optimset('TolX', 1e-15);
            for i = 1:obj.pH
                rn(i) = fzero(@(beta_)ZHn(obj, beta_), interval(i), options);
            end
            beta_ = rn(diff(rn)>1e-10);        % only keep unique roots
            beta_ = beta_(beta_ > 1e-10); % remove zeros                   
        end 
        function beta_ = computeBetaFH(obj)
            interval = linspace(0, obj.pfFH, obj.pFH_);   % interval/spacing 
                                                     % of root calculation
            rm = NaN*ones(obj.pFH_, 1);              % roots initialization
            options = optimset('TolX', 1e-15);
            for i = 1:obj.pFH_
                rm(i) = fzero(@(beta_) ZHn(obj, beta_), interval(i), options);
            end
            beta_ = rm(diff(rm) > 1e-10);     % only keep unique roots
            beta_ = beta_(beta_ > 1e-10);   % remove zeros   
        end 
        % computation of eta values
        function eta_ = computeEtaH(obj)
            interval = linspace(0, obj.qfH, obj.qH);    % interval/spacing 
                                                        % of root calculation
            rm = NaN*ones(obj.qH, 1);                   % roots initialization
            options = optimset('TolX', 1e-15);
            for i = 1:obj.qH
                    rm(i) = fzero(@(eta_) ZHm(obj, eta_), ...
                                                    interval(i), options);
            end
            eta_ = rm(diff(rm) > 1e-10);       % only keep unique roots
            eta_(1) = 0;
%             obj.etaHStatic = eta_;   
        end    
        % computation of fourier coefficients
        function [Cn_, beta_] = computeFHhatCn(obj, beta_)
            % compute beta and eta values if not already populated
            if length(beta_) < obj.niFH
                ni_ = length(beta_);
            else
                ni_ = obj.niFH;
            end
            CnTemp = NaN*ones(1, ni_);
            betaTemp = NaN*ones(1, ni_);            
            for j = 1:ni_
                CnTemp(j) = FourierCoefficientFH1(obj, beta_(j));
                betaTemp(j) = beta_(j);
            end
            cGet_ = abs(CnTemp) > obj.climFH;
            Cn_ = CnTemp(cGet_);
            beta_ = betaTemp(cGet_);                                     
        end
        function [Cnm_, beta_, eta_] = computeCnmH(obj, IC, FH, z_, r_, beta_, eta_)
            % compute beta and eta values if not already populated
%             if isempty(obj.etaHStatic)
%                 computeEtaH(obj);
%             end
            if length(beta_) < obj.niH 
                ni_ = length(beta_);
            else
                ni_ = obj.niH;
            end
            if length(eta_) < obj.miH
                mi_ = length(eta_);
            else
                mi_ = obj.miH;
            end
            CnmTemp = NaN*ones(ni_, mi_);
            betaTemp = NaN*ones(ni_, mi_);
            etaTemp = NaN*ones(ni_, mi_);            
            for i = 1:ni_
                for j = 1:mi_
                    CnmTemp(i, j) = FourierCoefficientThetaHhat(obj, ...
                        IC, FH, z_, r_, beta_(i), eta_(j));
                    betaTemp(i, j) = beta_(i);
                    etaTemp(i, j) = eta_(j);
                end
            end
            cGet_ = abs(CnmTemp) > obj.climH;
            Cnm_ = CnmTemp(cGet_);
            beta_ = betaTemp(cGet_);
            eta_ = etaTemp(cGet_); 
%             w = obj.betaH*pi/max(obj.betaH); 
%             sigmaN = sin(w)./w;
            sigmaN = ones(length(beta_), 1);
            Cnm_ = sigmaN.*Cnm_;              
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % conduction temperature solution for discharge mode
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function thetaS_ = iterateThetaS(obj, thetaC_, zcenter_, thetaT_, rtop_, thetaI_, IC_, z_, r_)      
            % compute auxilary problem solution 
            % compute homogeneous filters
            FD_ = computeFD(obj, z_, r_);
            % eigenvalues and Fourier coeffiecients
            [Cnm_, beta_, eta_] = computeGCnm(obj, IC_, FD_, thetaI_, z_, r_);
            % initial condition contribution
            thetaIC_ = computeThetaIC(obj, thetaI_, z_, r_, Cnm_, beta_, eta_);
            % compute boundary condition contributions
            [thetaBC1_, thetaBC3_] = computeThetaBCNow(obj, thetaC_, zcenter_, thetaT_, rtop_, thetaI_, z_, r_, beta_, eta_);            
            % sum individual solution terms
            thetaS_ = thetaIC_ + thetaBC1_ + thetaBC3_ - FD_;
        end          
        % z-dimension temperature BVP
        function xi = Xn(~, zbar_, beta_)
            % eigenfunction for z-dimension BVP
            xi = cos(beta_*zbar_);
        end
        function ni = Nn(obj, beta_)
            % norm for z-dimension BVP
            ni = 0.5*(obj.ztop*(beta_.^2 + obj.Bi1^2) + obj.Bi1)./ ...
                                                (beta_.^2 + obj.Bi1^2);
        end
        function zi = Zn(obj, beta_)
            % transendental function for beta values
            zi = beta_.*sin(obj.ztop*beta_) - obj.Bi1*cos(obj.ztop*beta_);
        end     
        function xi = XFn(obj, zbar_, eta_)
            % z-BVP eigenfunction for filter 2
            xi = eta_.*(sinh(eta_*obj.ztop).*sinh(eta_*zbar_) - ...
                        cosh(eta_*obj.ztop).*cosh(eta_*zbar_)) + ...
               obj.Bi1*(cosh(eta_*obj.ztop).*sinh(eta_*zbar_) - ...
                        sinh(eta_*obj.ztop).*cosh(eta_*zbar_));
        end
        % r-dimension temperature BVP
        function si = S(obj, eta_)
            si = -eta_.*bessely(1, eta_*obj.b) + ...
                obj.Bi4*bessely(0, eta_*obj.b);
        end
        function ui = U(obj, eta_)
            ui = -eta_.*besselj(1, eta_*obj.a0) - ...
                obj.Bi3*besselj(0, eta_*obj.a0);
        end
        function vi = V(obj, eta_)
            vi = -eta_.*besselj(1, eta_*obj.b) + ...
                obj.Bi4*besselj(0, eta_*obj.b);
        end
        function wi = W(obj, eta_)
            wi = -eta_.*bessely(1, eta_*obj.a0) - ...
                obj.Bi3*bessely(0, eta_*obj.a0);
        end
        function b1i = B1(obj, eta_)
            b1i = obj.Bi3^2 + eta_.^2;
        end
        function b2i = B2(obj, eta_) 
            b2i = obj.Bi4^2 + eta_.^2;
        end
        function xi = Xm(obj, rbar_, eta_)
            % eigenfunction for r-dimension BVP
            xi = besselj(1, eta_*obj.b)*bessely(0, eta_*rbar_) - ...
                 besselj(0, eta_*rbar_)*bessely(1, eta_*obj.b);
        end
        function ni = Nm(obj, eta_)
            % norm for r-dimension BVP
            ni = (2/pi^2)*(U(obj, eta_).^2 - B1(obj, eta_).* ...
                   besselj(1, eta_*obj.b).^2)./(eta_.^2.*U(obj, eta_).^2);
        end
        function zi = Zm(obj, eta_)
            % transendental function for eta values
            zi = W(obj, eta_).*besselj(1, eta_*obj.b) - ...
                 U(obj, eta_).*bessely(1, eta_*obj.b);
        end   
        function [xi, xip] = XFm(obj, rbar_, beta_)
            % eigenfunction for r-component of filter 1
            xi = obj.Bi3*(besselk(0, beta_*rbar_).* ...
                          besseli(0, beta_*obj.a0) - ...
                          besselk(0, beta_*obj.a0).* ...
                          besseli(0, beta_*rbar_)) + ...
                 beta_*(besselk(0, beta_*rbar_).* ...
                        besseli(1, beta_*obj.a0) - ...
                        besselk(1, beta_*obj.a0).* ...
                        besseli(0, beta_*rbar_));
            xip = -beta_*(obj.Bi3*(besselk(1, beta_*rbar_).* ...
                                   besseli(0, beta_*obj.a0) + ...
                                   besselk(0, beta_*obj.a0).* ...
                                   besseli(1, beta_*rbar_)) + ...
                          beta_*(besselk(1, beta_*rbar_).* ...
                                 besseli(1, beta_*obj.a0) + ...
                                 besselk(1, beta_*obj.a0).* ...
                                 besseli(1, beta_*rbar_)));            
        end
        % combined homogeneous temperature solution           
        function c = FourierCoefficientG(obj, IC_, FD_, thetaI_, z_, r_, beta_, eta_, R_)
            % Fourier coefficients
            C_num = computeICIntegral(obj, IC_, FD_, thetaI_, z_, r_, beta_, R_);
            c = C_num/(Nn(obj, beta_)*Nm(obj, eta_));
        end   
        function c = FourierCoefficientFD1(obj, beta_)
            % Fourier coefficients
            [~, Rpb] = XFm(obj, obj.b, beta_);
            C_num = mean(mean(obj.g4))*sin(beta_*obj.ztop)/(beta_*Rpb);
            c = C_num/Nn(obj, beta_);  
        end
        function c = FourierCoefficientFD2(obj, eta_)
            % Fourier coefficients
            C_num = obj.g2*(besselj(1, eta_*obj.b).* ...
                                   (obj.b*bessely(1, eta_*obj.b) - ...
                                    obj.a0*bessely(1, eta_*obj.a0)) - ...
                            bessely(1, eta_*obj.b).* ...
                                   (obj.b*besselj(1, eta_*obj.b) - ...
                                    obj.a0*besselj(1, eta_*obj.a0)))/ ...
                    (eta_.^2.*(eta_.*sinh(eta_*obj.ztop) + ...
                               obj.Bi1*cosh(eta_*obj.ztop)));
            c = C_num/Nm(obj, eta_);  
        end
        function thetaIC_ = computeThetaIC(obj, thetaI_, z_, r_, Cnm_, beta_, eta_)
            % compute eigenvalues and coefficients if not yet populated         
            n = length(z_); m = length(r_);
            % compute for current time step
            thetaIC_ = zeros(n, m);
            for i = 1:length(Cnm_)
                C = Cnm_(i);
                F = Xt(obj, obj.df, beta_(i), eta_(i));
                Z = Xn(obj, z_, beta_(i));  
                R = Xm(obj, r_, eta_(i));
                thetaIC_ = thetaIC_ + C*F*Z'*R;
            end 
            thetaIC_ = thetaIC_ + thetaI_;
        end
        function x = computeICIntegral(obj, IC_, FD_, thetaI_, z_, r_, beta_, R_)
            % computes integral component of the initial condition
            % contribution in Green's formula
            fr = simpsonIntegrator(obj, r_);
            fz = simpsonIntegrator(obj, z_);
            Ir = (R_.*r_.*(IC_ + FD_ - thetaI_))*fr';
            x = (Xn(obj, z_, beta_).*Ir')*fz';
        end          
        function FD_ = computeFD(obj, z_, r_)
            % computes the filtering function for the bin holding solution
            if isempty(obj.g2), initializeBaseSys(obj); end
            if isempty(obj.g4), initializeWallSys(obj); end
            beta_ = computeBetaFD(obj);
            eta_ = computeEtaFD(obj);
            [Cn_, beta_] = computeFDCn(obj, beta_);  
            [Cm_, eta_] = computeFDCm(obj, eta_);
            % compute FHhat
            n = length(z_); m = length(r_);
            FD1 = zeros(n, m); FD2 = zeros(n, m);
            for i = 1:length(Cn_)
                Z = Xn(obj, z_, beta_(i));
                R = XFm(obj, r_, beta_(i));
                FD1 = FD1 + Cn_(i)*Z'*R;               
            end
            for i = 1:length(Cm_)
                Z = XFn(obj, z_, eta_(i));
                R = Xm(obj, r_, eta_(i));
                FD2 = FD2 + Cm_(i)*Z'*R;               
            end
            FD_ = FD1 + FD2;
        end
        function [thetaBC1_, thetaBC3_] = computeThetaBCNow(obj, thetaC_, zcenter_, thetaT_, rtop_, thetaI_, z_, r_, beta_, eta_)
            % computes boundary integration value at current time step
            n = length(z_); m = length(r_);
            thetaBC1_ = zeros(n, m); thetaBC3_ = zeros(n, m);
            yt = obj.T2S(thetaT_, rtop_, r_, 'top'); 
            yc = obj.C2S(thetaC_, zcenter_, z_, 'top');
            fr = simpsonIntegrator(obj, r_);
            fz = simpsonIntegrator(obj, z_);
            % compute top and center boundary temperatures 
            g1_ = obj.Bi1*(yt - thetaI_);
            g3_ = obj.Bi3*(yc - thetaI_);
            if isempty(obj.g1)
                obj.g1 = obj.Bi1*(ones(size(g1_)) - thetaI_); 
            end
            if isempty(obj.g3)
                obj.g3 = obj.Bi3*(ones(size(g3_)) - thetaI_); 
            end
            for i = 1:length(beta_)
                RD_ = Xm(obj, r_, eta_(i));
                RC3_ = Xm(obj, obj.a0, eta_(i));
                Z = Xn(obj, z_, beta_(i))';
                C1 = Xn(obj, z_(end), beta_(i)) ...
                    /(Nn(obj, beta_(i))*Nm(obj, eta_(i)));
                C3 = obj.a0*RC3_ ...
                    /(Nn(obj, beta_(i))*Nm(obj, eta_(i)));
                IBC1_ = (obj.g1.*RD_.*r_)*fr';
                IBC3_ = (obj.g3(1:n).*Xn(obj, z_, beta_(i)))*fz';
                IBC1 = (g1_.*RD_.*r_)*fr';
                IBC3 = (g3_.*Xn(obj, z_, beta_(i)))*fz';
                p1 = (IBC1 - IBC1_)/obj.df; q1 = IBC1_;
                p3 = (IBC3 - IBC3_)/obj.df; q3 = IBC3_;
                lambda = sqrt(beta_(i)^2 + eta_(i)^2);
                F1 = q1*(1 - exp(-obj.df*lambda^2))/lambda^2 + ...
                     p1*(1 - (lambda^2*obj.df + 1)* ...
                      exp(-obj.df*lambda^2))/lambda^4;
                F3 = q3*(1 - exp(-obj.df*lambda^2))/lambda^2 + ...
                     p3*(1 - (lambda^2*obj.df + 1)* ...
                      exp(-obj.df*lambda^2))/lambda^4;
                ep1 = C1*F1*Z*RD_; 
                ep3 = C3*F3*Z*RD_; 
                thetaBC1_ = thetaBC1_ + ep1;
                thetaBC3_ = thetaBC3_ + ep3;
            end
            obj.g1 = g1_; obj.g3 = g3_;
        end
        % computation of beta values
        function beta_ = computeBeta(obj)
            interval = linspace(0, obj.pf, obj.p);   % interval/spacing 
                                                     % of root calculation
            rn = NaN*ones(obj.p, 1);                 % roots initialization
            options = optimset('TolX', 1e-15);
            for i = 1:obj.p
                rn(i) = fzero(@(beta_)Zn(obj, beta_), interval(i), options);
            end
            beta_ = rn(diff(rn) > 1e-10);         % only keep unique roots
            beta_ = beta_(beta_ > 1e-10);         % remove zeros                              
        end       
        function beta_ = computeBetaFD(obj)
            interval = linspace(0, obj.pfFD, obj.pFD);   % interval/spacing 
                                                     % of root calculation
            rn = NaN*ones(obj.pFD, 1);                 % roots initialization
            options = optimset('TolX', 1e-15);
            for i = 1:obj.pFD
                rn(i) = fzero(@(beta_)Zn(obj, beta_), interval(i), options);
            end
            beta_ = rn(diff(rn)>1e-10);         % only keep unique roots
            beta_ = beta_(beta_ > 1e-10); % remove zeros                              
        end   
        % computation of eta values
        function eta_ = computeEta(obj)
            interval = linspace(1, obj.qf, obj.q);   % interval/spacing 
                                                     % of root calculation
            rm = NaN*ones(obj.q, 1);                 % roots initialization
            options = optimset('TolX', 1e-15);
            for i = 1:obj.q
                rm(i) = fzero(@(eta_) Zm(obj, eta_), interval(i), options);
            end
            eta_ = rm(diff(rm) > 1e-10);         % only keep unique roots
            eta_ = eta_(eta_ > 1e-10);           % remove zeros  
%             % save eta and eta dependencies for future reference
%             obj.etaStatic = obj.eta;
%             obj.RStatic = cell(length(obj.eta), 1);
%             obj.RC3Static = NaN*ones(length(obj.eta), 1);
%             obj.RC4Static = NaN*ones(length(obj.eta), 1);
%             for i = 1:length(obj.eta)            
%                 obj.RStatic{i} = Xm(obj, obj.rbar, obj.eta(i));
%                 obj.RC3Static(i) = Xm(obj, obj.a0, obj.eta(i));
%                 obj.RC4Static(i) = Xm(obj, obj.b, obj.eta(i));
%             end   
        end       
        function eta_ = computeEtaFD(obj)
            interval = linspace(1, obj.qfFD, obj.qFD);   % interval/spacing 
                                                     % of root calculation
            rm = NaN*ones(obj.qFD, 1);                 % roots initialization
            options = optimset('TolX', 1e-15);
            for i = 1:obj.qFD
                rm(i) = fzero(@(eta_) Zm(obj, eta_), interval(i), options);
            end
            eta_ = rm(diff(rm)>1e-10);          % only keep unique roots
            eta_ = eta_(eta_ > 1e-10);          % remove zeros  
            % save eta and eta dependencies for future reference
%             obj.etaFDStatic = obj.etaFD;
%             obj.RFDStatic = cell(length(obj.etaFD), 1);
%             for i = 1:length(obj.etaFD)            
%                 obj.RFDStatic{i} = Xm(obj, obj.rbar, obj.etaFD(i));
%             end 
        end   
        % computation of fourier coefficients
        function [Cnm_, beta_, eta_] = computeGCnm(obj, IC_, FD_, thetaI_, z_, r_)
            % compute new beta values
            beta_ = computeBeta(obj);
            % compute eta values if not already populated
            eta_ = computeEta(obj);
            if length(beta_) < obj.ni
                obj.ni = length(beta_);
            end
            if length(eta_) < obj.mi
                obj.mi = length(eta_);
            end
            CnmTemp = NaN*ones(obj.ni, obj.mi);
            betaTemp = NaN*ones(obj.ni, obj.mi);
            etaTemp = NaN*ones(obj.ni, obj.mi);            
%             RC3Temp = NaN*ones(obj.ni, obj.mi);
%             RC4Temp = NaN*ones(obj.ni, obj.mi);
%             RTemp = cell(obj.ni, obj.mi);
            for i = 1:obj.ni
                for j = 1:obj.mi
                    R_ = Xm(obj, r_, eta_(j));
                    CnmTemp(i, j) = FourierCoefficientG(obj, IC_, FD_, ...
                        thetaI_, z_, r_, beta_(i), eta_(j), R_);
                    betaTemp(i, j) = beta_(i);
                    etaTemp(i, j) = eta_(j);
%                     RTemp{i, j} = obj.RStatic{j};
%                     RC3Temp(i, j) = obj.RC3Static(j);
%                     RC4Temp(i, j) = obj.RC4Static(j);
                end
            end
            cGet_ = abs(CnmTemp) > obj.clim;
            Cnm_ = CnmTemp(cGet_);
            beta_ = betaTemp(cGet_);
%             w = obj.beta*pi/max(obj.beta); 
%             sigmaN = sin(w)./w;
            sigmaN = ones(length(beta_), 1);
            Cnm_ = sigmaN.*Cnm_;
            eta_ = etaTemp(cGet_);
%             obj.RD = RTemp(obj.cGet);
%             obj.RC3 = RC3Temp(obj.cGet);
%             obj.RC4 = RC4Temp(obj.cGet);               
        end
        function [Cn_, beta_] = computeFDCn(obj, beta_)
            % compute beta and eta values if not already populated
            if length(beta_) < obj.niFD
                obj.niFD = length(beta_);
            end
            CnTemp = NaN*ones(1, obj.niFD);
            betaTemp = NaN*ones(1, obj.niFD);            
            for j = 1:obj.niFD
                CnTemp(j) = FourierCoefficientFD1(obj, beta_(j));
                betaTemp(j) = beta_(j);
            end
            cGet_ = abs(CnTemp) > obj.climFD;
            Cn_ = CnTemp(cGet_);
            beta_ = betaTemp(cGet_);               
        end
        function [Cm_, eta_] = computeFDCm(obj, eta_)
            % compute beta and eta values if not already populated
            if length(eta_) < obj.miFD
                obj.miFD = length(eta_);
            end
            CmTemp = NaN*ones(1, obj.miFD);           
            for j = 1:obj.miFD
                CmTemp(j) = FourierCoefficientFD2(obj, eta_(j));
            end
            cGet_ = abs(CmTemp) > obj.climFD;
            Cm_ = CmTemp(cGet_);
            eta_ = eta_(cGet_);
%             obj.RFD = obj.RFDStatic(cGet_)                                   
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % top boundary velocity solution for discharge mode
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function ubar_ = computeUbar(obj, rtop_)
            % computes velocity profile on top boundary
            ubar_ = -obj.Qc./(2*pi*rtop_*obj.h);
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
            Ls = diag(ones(length(ubar_)-1, 1), -1);
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
            wbar_ = -obj.Qc/(pi*(obj.a0)^2)*ones(size(zcenter_));             
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
            OmegaC3_ = 1/obj.dzc^2 + obj.Pe ...
                   *repmat(wbar_', m, 1)/(2*obj.dzc);
            OmegaC4_ = 1/obj.drhat^2 ...
                        - 1/(2*repmat(rhat_', n, 1)*obj.drhat);
            OmegaC5_ = 1/obj.dzc^2 - obj.Pe ...
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
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % storage bin quiescent air temperature
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function computeThetaA(obj, thetaP, z_, r_)
            % computes the air temperature at the current iteration
            t = obj.Fo2t(obj.df, 1); n = length(z_); m = length(r_);
            % heat loss (W) at the particle boundary
            if isempty(obj.qTopP)
                obj.qTopP = cell(2, 1);
                obj.qTopP{1} = zeros(1, m);
            else
                obj.qTopP{1} = obj.qTopP{2};
            end
            fr = simpsonIntegrator(obj, r_);
            obj.qTopP{2} = -obj.kp*(thetaP(end, :) - thetaP(end-1, :))./ ...
                         (z_(end) - z_(end-1))*(obj.T0 - obj.Tinf)/obj.Hp;
            qTopP_ = 2*pi*fr*(obj.qTopP{1}.*r_)'*obj.Hp^2;
            qTopP__ = 2*pi*fr*((obj.qTopP{2} - obj.qTopP{1}) ...
                                                     ./t.*r_)'*obj.Hp^2; 
            % heat loss (W) at the top bin boundary
            At = pi*obj.bp^2*obj.Hp^2;
            qTop_ = obj.qLossT(1)*At;
            qTop__ = (obj.qLossT(2) - obj.qLossT(1))/t*At;
            % heat loss (W) at wall boundary
            [~, zwi] = min(abs(obj.ztop - obj.zbarW));
            zw = obj.zbarW(zwi:end);
            fz = simpsonIntegrator(obj, zw);
            qWall_ = 2*pi*obj.bp*fz*obj.qWall{1, 1}(zwi:end)*obj.Hp^2;
            qWall__ = 2*pi*obj.bp*fz*(obj.qWall{2, 1}(zwi:end) - ...
                                   obj.qWall{1, 1}(zwi:end))/t*obj.Hp^2;
            % compute thetaA
            cp_ = 1005;     % (J/kgK)
            rho_ = 1.225;   % (m3/kg)
            V = pi*obj.b^2*(1 - obj.ztop)*obj.Hp^3;
            C_ = cp_*rho_*V;            
            TA = obj.theta2T(obj.thetaA) + t*(qTopP_ - qTop_ - qWall_)/C_ ...
                + 0.5*t^2*(qTopP__ - qTop__ - qWall__)/C_;
            obj.thetaA = obj.T2theta(TA);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % storage bin top, base and wall RC model
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function computeUtop(obj)
            % computes the overall heat transfer coefficient for the base
            % using the insulation cell array
            R = 0;
            for i = 1:size(obj.roofInsulation, 1)
                t = abs(obj.roofInsulation{i, 2}(2) - ...
                                            obj.roofInsulation{i, 2}(1));
                k_ = obj.roofInsulation{i, 3};
                R = R + t/(pi*(obj.bp*obj.Hp)^2*k_);
            end
            R = R + 1/(pi*(obj.bp*obj.Hp)^2*obj.hInf);            
            obj.hpTop = 1/(pi*(obj.bp*obj.Hp)^2*R);
        end
        function computeUbase(obj)
            % computes the overall heat transfer coefficient for the base
            % using the insulation cell array
            R = 0;
            for i = 1:size(obj.baseInsulation, 1)
                t = abs(obj.baseInsulation{i, 2}(2) - ...
                                            obj.baseInsulation{i, 2}(1));
                k_ = obj.baseInsulation{i, 3};
                R = R + t/(pi*(obj.bp*obj.Hp)^2*k_);
            end
            R = R + 1/(pi*(obj.bp*obj.Hp)^2*obj.hInf);            
            obj.hp2 = 1/(pi*(obj.bp*obj.Hp)^2*R);
        end
        function computeUwall(obj)
            % computes the overall heat transfer coefficient for the wall
            % using the insulation cell array
            R = 0;
            for i = 1:size(obj.wallInsulation, 1)
                r1 = obj.wallInsulation{i, 2}(1);
                r2 = obj.wallInsulation{i, 2}(2);
                k_ = obj.wallInsulation{i, 3};
                R = R + log(r2/r1)/(2*pi*obj.Hp*k_);
            end
            R = R + 1/(2*pi*r2*obj.Hp*obj.hInf);
            obj.hp4 = 1/(2*pi*obj.bp*obj.Hp^2*R);
        end
        function initializeTopSys(obj)
            % computes the linear RC system for the storage tank base
            N = size(obj.roofInsulation, 1); 
            obj.Rtop = {}; obj.Ctop = {}; 
            % insulation layer capacitance and resistance
            for i = 1:N
                obj.Rtop{i, 1} = obj.roofInsulation{i, 1};
                obj.Ctop{i, 1} = obj.roofInsulation{i, 1};
                t = obj.roofInsulation{i, 2}(2) - ...
                                          obj.roofInsulation{i, 2}(1);
                k_ = obj.roofInsulation{i, 3};
                rho_ = obj.roofInsulation{i, 4};
                c_ = obj.roofInsulation{i, 5};
                obj.Rtop{i, 2} = t/(pi*(obj.bp*obj.Hp)^2*k_);
                obj.Ctop{i, 2} = rho_*c_*t;
            end
            % convective resistance
            obj.Rtop{N + 1, 1} = 'convection';
            obj.Rtop{N + 1, 2} = 1/(pi*(obj.bp*obj.Hp)^2*obj.hInf);
            % time constants
            obj.tauTop = NaN*ones(N, 2);
            obj.tauTop(:, 1) = obj.Hp^2./(obj.alphapPacked* ...
                              [obj.Ctop{:, 2}].*[obj.Rtop{1:end-1, 2}]);
            obj.tauTop(:, 2) = obj.Hp^2./(obj.alphapPacked* ...
                                [obj.Ctop{:, 2}].*[obj.Rtop{2:end, 2}]);
            % state-space system
            obj.Atop = spdiags([-(obj.tauTop(:, 1) + obj.tauTop(:, 2)), ...
                        obj.tauTop(:, 1), obj.tauTop(:, 2)], ...
                        [0, 1, -1], N, N)';
            obj.Btop = zeros(N, 1); obj.Btop(1) = obj.tauTop(1, 1);
            obj.Ctop = [zeros(1, N); eye(N)]; 
            obj.Ctop(1) = -(obj.T0 - obj.Tinf)/obj.Rtop{1, 2};
            obj.Dtop = zeros(N+1, 1);
            obj.Dtop(1) = (obj.T0 - obj.Tinf)/obj.Rtop{1, 2};
            obj.topSys = ss(full(obj.Atop), obj.Btop, obj.Ctop, obj.Dtop); 
            % initialize state variables and boundary condition
            if isempty(obj.thetaTop), obj.thetaTop = zeros(N, 1); end
            computeTopSys(obj);
        end
        function initializeBaseSys(obj)
            % computes the linear RC system for the storage tank base
            N = size(obj.baseInsulation, 1); 
            obj.Rbase = cell(N+1, 2); obj.CapBase = cell(N, 2); 
            % insulation layer capacitance and resistance
            for i = 1:N
                obj.Rbase{i, 1} = obj.baseInsulation{i, 1};
                obj.CapBase{i, 1} = obj.baseInsulation{i, 1};
                t = abs(obj.baseInsulation{i, 2}(2) - ...
                                 obj.baseInsulation{i, 2}(1));
                k_ = obj.baseInsulation{i, 3};
                rho_ = obj.baseInsulation{i, 4};
                c_ = obj.baseInsulation{i, 5};
                obj.Rbase{i, 2} = t/(pi*(obj.b*obj.Hp)^2*k_);
                obj.CapBase{i, 2} = rho_*c_*t*(obj.b*obj.Hp)^2;
            end
            % convective resistance
            obj.Rbase{N + 1, 1} = 'convection';
            obj.Rbase{N + 1, 2} = 1/(pi*(obj.b*obj.Hp)^2*obj.hInf);
            % time constants
            obj.tauBase = NaN*ones(N, 2);
            obj.tauBase(:, 1) = obj.Hp^2./(obj.alphapPacked* ...
                              [obj.CapBase{:, 2}].*[obj.Rbase{1:end-1, 2}]);
            obj.tauBase(:, 2) = obj.Hp^2./(obj.alphapPacked* ...
                                [obj.CapBase{:, 2}].*[obj.Rbase{2:end, 2}]);
            % state-space system
            obj.Abase = spdiags([-(obj.tauBase(:, 1) + obj.tauBase(:, 2)), ...
                        obj.tauBase(:, 1), obj.tauBase(:, 2)], ...
                        [0, 1, -1], N, N)';
            obj.Bbase = zeros(N, 1); obj.Bbase(1) = obj.tauBase(1, 1);
            obj.Cbase = [zeros(1, N); eye(N)]; 
            obj.Cbase(1) = -(obj.T0 - obj.Tinf)/(obj.Rbase{1, 2}*pi*(obj.b*obj.Hp)^2);
            obj.Dbase = zeros(N+1, 1);
            obj.Dbase(1) = (obj.T0 - obj.Tinf)/(obj.Rbase{1, 2}*pi*(obj.b*obj.Hp)^2); 
            % initialize state variables and boundary condition
            if isempty(obj.thetaBase), obj.thetaBase = zeros(N, 1); end
        end
        function initializeWallSys(obj)
            % computes the linear RC system for the storage tank base
            N = size(obj.wallInsulation, 1); 
            obj.Rwall = cell(N+1, 2); obj.CapWall = cell(N, 2); 
            % insulation layer capacitance and resistance
            for i = 1:N
                obj.Rwall{i, 1} = obj.wallInsulation{i, 1};
                obj.CapWall{i, 1} = obj.wallInsulation{i, 1};
                r1 = obj.wallInsulation{i, 2}(1);
                r2 = obj.wallInsulation{i, 2}(2);
                k_ = obj.wallInsulation{i, 3};
                rho_ = obj.wallInsulation{i, 4};
                c_ = obj.wallInsulation{i, 5};
                obj.Rwall{i, 2} = log(r2/r1)/(2*pi*obj.Hp*k_);
                obj.CapWall{i, 2} = rho_*c_*obj.Hp*pi*(r2^2 - r1^2);
            end
            % convective resistance
            obj.Rwall{N + 1, 1} = 'convection';
            obj.Rwall{N + 1, 2} = 1/(2*pi*r2*obj.Hp*obj.hInf);
            % time constants
            obj.tauWall = NaN*ones(N, 2);
            obj.tauWall(:, 1) = obj.Hp^2./(obj.alphapPacked* ...
                                [obj.CapWall{:, 2}].*[obj.Rwall{1:end-1, 2}]);
            obj.tauWall(:, 2) = obj.Hp^2./(obj.alphapPacked* ...
                                [obj.CapWall{:, 2}].*[obj.Rwall{2:end, 2}]);
            % state-space system
            obj.Awall = spdiags([-(obj.tauWall(:, 1) + obj.tauWall(:, 2)), ...
                        obj.tauWall(:, 1), obj.tauWall(:, 2)], ...
                        [0, 1, -1], N, N)';
            obj.Bwall = zeros(N, 1); obj.Bwall(1) = obj.tauWall(1, 1);
            obj.Cwall = [zeros(1, N); eye(N)];  
            for i = 1:N
                obj.Cwall(i, i) = -(obj.T0 - obj.Tinf)/(obj.Rwall{i, 2}*2*pi*obj.wallInsulation{i, 2}(1)*obj.Hp);
                if i ~= 1
                    obj.Cwall(i, i-1) = (obj.T0 - obj.Tinf)/(obj.Rwall{i, 2}*2*pi*obj.wallInsulation{i, 2}(1)*obj.Hp);
                end
            end
            obj.Dwall = zeros(2*N, 1);
            obj.Dwall(1) = (obj.T0 - obj.Tinf)/(obj.Rwall{1, 2}*2*pi*obj.b*obj.Hp^2);
            % initialize state variables and boundary condition   
            if isempty(obj.thetaWall), obj.thetaWall = zeros(N, 1); end
        end
        function y = computeTopSys(obj, u, Fo_)
            % computes the heat flux leaving the
            if nargin < 2, u = zeros(2, 1); end
            if nargin < 3, Fo_ = linspace(0, obj.df, length(u)); end
            if isempty(obj.qLossT), obj.qLossT = zeros(2, 1); end
            y = lsim(obj.topSys, u, Fo_, obj.thetaTop);
            obj.qLossT(1) = obj.qLossT(2); obj.qLossT(2) = y(end, 1); 
            obj.thetaTop = y(end, 2:end);
        end
        function y = computeBaseSys(obj, u, Fo_)
            % computes the heat flux leaving the
            if nargin < 2, u = 1; end
            if nargin < 3, Fo_ = obj.df; end
            N = size(obj.baseInsulation, 1);
            b_ = zeros(N, 1);
            b_(1) = obj.tauBase(1, 1)*u;
            y_ = expm(Fo_.*[full(obj.Abase) eye(N); zeros(N) zeros(N)]) ...
                *[obj.thetaBase; b_];
            y = y_(1:N);
            obj.thetaBase = y;
            obj.qLossB = (obj.T0 - obj.Tinf)/(obj.Rbase{1, 2}*pi*(obj.b*obj.Hp)^2)*(u - y(1));
            obj.g2 = -obj.qLossB*obj.H_/(obj.k*(obj.T0 - obj.Tinf));           
        end
        function y = computeWallSys(obj, u, Fo_)
            % computes the heat flux leaving the wall
            if nargin < 2, u = 1; end
            if nargin < 3, Fo_ = obj.df; end
            N = size(obj.wallInsulation, 1);
            b_ = zeros(N, 1); 
            b_(1) = obj.tauWall(1, 1)*u;
            y_ = expm(Fo_.*[full(obj.Awall) eye(N); zeros(N) zeros(N)]) ...
                *[obj.thetaWall; b_];
            y = y_(1:N);
            obj.thetaWall = y;
            obj.qLossW = (obj.T0 - obj.Tinf)/(obj.Rwall{1, 2}*2*pi*obj.b*obj.Hp^2)*(u - y(1));
            obj.g4 = obj.qLossW*obj.H_/(obj.k*(obj.T0 - obj.Tinf));
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Domain modification or data filling
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function updateDomain(obj)
            % updates domain according to mass accounting
            storeTotalMass(obj);
            obj.HNow = obj.HNow - obj.dH;
            obj.ztop = obj.HNow/obj.H_;
            fitzMesh(obj);
            fitTheta(obj);
        end
        function updateDomainC(obj)
            % updates domain according to mass accounting
            storeTotalMass(obj);
            obj.HNow = obj.HNow + obj.dHCh;
            obj.ztop = obj.HNow/obj.H_;
            fitzMeshC(obj);
            fitThetaC(obj);
        end
        function fitzMeshC(obj)
            % fits z meshes to be between 0 and ztop, resizing dzbar and
            % dzc accordingly. Adjusts the length according to the time
            % step relative to the set update frequency, modZ.
            [~, n] = min(abs(obj.FoNow - obj.Fo));
            if mod(n, obj.modZH) == 0
                obj.nzH = length(obj.zbarH) + 1; 
            else
                obj.nzH = length(obj.zbarH);
            end
            [obj.zbarH, obj.dzH] = ...
                              nodeGen(obj, [0, obj.ztop], obj.nzH);         
        end
        function fitzMesh(obj)
            % fits z meshes to be between 0 and ztop, resizing dzbar and
            % dzc accordingly. Adjusts the length according to the time
            % step relative to the set update frequency, modZ.
%             [~, n] = min(abs(obj.FoNow - obj.Fo));
%             if mod(n, obj.modZS) == 0
%                 obj.nzbar = length(obj.zbar) - 1; 
%             else
%                 obj.nzbar = length(obj.zbar);
%             end
%             if mod(n, obj.modZC) == 0
%                 nzc = length(obj.zcenter) - 1;
%             else
%                 nzc = length(obj.zcenter);
%             end
            [obj.zbar, obj.dzbar] = ...
                            nodeGen(obj, [0, obj.ztop], obj.nzbar);
            obj.zcenter = linspace(0, obj.ztop, nzc);
            obj.dzc = obj.zcenter(2);
            computeWbar(obj);
            computeUbar(obj);            
        end
        function fitThetaC(obj)
            % fits/updates theta matrices to match current mesh
            move = abs(size(obj.thetaH, 1) - length(obj.zbarH));
            if move > 0             
                % update temperature 
                for n = 1:move
                    obj.thetaH = [obj.thetaH; ...
                           obj.thetaCi*ones(1, length(obj.thetaH(1, :)))];
                    % compute temperature distribution for mixing depth
                    setMixingProfileC(obj);
                    % set new initial condition for next time step
                    obj.rhoH = obj.thetaH; 
                end
            end         
        end
        function setMixingProfileC(obj)
            % computes bulk temperature/s for the mixing depth for the 
            % current time step in charging mode
            [~, im] = min(abs(obj.zbarH - (obj.ztop - obj.deltaM)));
            zm = obj.zbarH(im:end);
            if length(zm) < 3
                im = im - 3;
                zm = obj.zbarH(im:end);
            end
            obj.thetaMC = mean(obj.thetaH(im:end, :));
            obj.thetaH(im:end, :) = repmat(obj.thetaMC, length(zm), 1);
        end
        function fitTheta(obj)
            % fits/updates theta matrices to match current mesh
            move = abs(size(obj.thetaS, 1) - length(obj.zbar));
            movec = abs(size(obj.thetaC, 1) - length(obj.zcenter));
            if move > 0             
                % update temperature 
                for n = 1:move
                    for i = 1:length(obj.zhat)-1
                        obj.thetaT(i, :) = obj.thetaT(i+1, :);                    
                    end
                    obj.thetaT(end, :) = obj.S2T;
                    obj.thetaS = obj.thetaS(1:end-1, :); 
                    obj.FS0 = obj.FS0(1:end-1, :);
                end
            end
            if movec > 0
                for n = 1:movec
                    obj.thetaC = obj.thetaC(1:end-1, :);
                    obj.scPe = obj.scPe(1:end-1, :);
                    obj.scQc = obj.scQc(1:end-1, :);
                    obj.sca = obj.sca(1:end-1, :);                    
                end
            end             
        end
        function theta_ = patchTheta(obj, thetaS_, thetaT_, thetaC_, z_, r_, zhat_, zcenter_, rtop_, rhat_)
            % combines thetaS, thetaT and thetaC for current time step
            ns = length(zhat_); nl = length(z_); n = ns + nl - 1;
            ms = length(rhat_); ml = length(r_); m = ms + ml - 1;
            theta_ = NaN*ones(n, m);
            % overlay current temperature in stagnant region
            theta_(1:nl, ms:m) = thetaS_;
            % overlay current temperature in top region
            thetaTMatch = obj.T2S(thetaT_, rtop_, r_, 'full');
            theta_(nl+1:n, ms:m) = flipud(thetaTMatch(1:end-1, :));
            % overlay current temperature in center channel
            thetaCMatch = obj.C2S(thetaC_, zcenter_, z_, 'full');
            theta_(1:nl, 1:ms-1) = flipud(thetaCMatch(:, 1:end-1));
            % overlay bulk junction temperature in top-center junction
            theta_(nl+1:n, 1:ms-1) = mean(obj.thetaChat); 
        end
        function patchThetaC(obj, k_, verify_conservation)
            % save thetaH array and domain info
            if mod(k_-1, obj.ls) == 0 && k_ ~= 1
                saveTheta(obj, k_-1);
                obj.theta = {};                
            end
            % normalize time step
            kn_ = mod(k_-1, obj.ls) + 1;                        
            % append to theta cell array that will be saved to drive
            obj.theta{kn_, 1} = obj.thetaH;
            % store corresponding mesh arrays
            obj.theta{kn_, 2} = obj.zbarH;
            obj.theta{kn_, 3} = obj.rbarH; 
            obj.theta{kn_, 4} = k_;
            obj.theta{kn_, 5} = obj.Fo(k_);
            obj.theta{kn_, 7} = obj.thetaW;
            obj.theta{kn_, 8} = obj.thetaBase;
            obj.theta{kn_, 9} = obj.qWall;
            obj.theta{kn_, 10} = obj.thetaA;
            obj.theta{kn_, 11} = obj.qLossW;
            obj.theta{kn_, 12} = obj.qLossB;
            obj.theta{kn_, 13} = obj.qLossT;
            obj.theta{kn_, 14} = obj.qTopP;
            if verify_conservation
                computeEnergyH(obj);
                obj.theta{kn_, 15} = obj.energy;
            end
            % save if last time step
            if k_ == length(obj.Fo)
                saveTheta(obj, k_);
            end
        end
        function patchThetaH(obj, k_, verify_conservation)
            % save thetaH array and domain info
            if mod(k_-1, obj.ls) == 0 && k_ ~= 1
                saveTheta(obj, k_-1);
                obj.theta = {};                
            end
            % normalize time step
            kn_ = mod(k_-1, obj.ls) + 1;                        
            % append to theta cell array that will be saved to drive
            obj.theta{kn_, 1} = obj.thetaH;
            % store corresponding mesh arrays
            obj.theta{kn_, 2} = obj.zbarH;
            obj.theta{kn_, 3} = obj.rbarH; 
            obj.theta{kn_, 4} = k_;
            obj.theta{kn_, 5} = obj.Fo(k_);
            obj.theta{kn_, 7} = obj.thetaW;
            obj.theta{kn_, 8} = obj.thetaBase;
            obj.theta{kn_, 9} = obj.qWall;
            obj.theta{kn_, 10} = obj.thetaA;
            obj.theta{kn_, 11} = obj.qLossW;
            obj.theta{kn_, 12} = obj.qLossB;
            obj.theta{kn_, 13} = obj.qLossT;
            obj.theta{kn_, 14} = obj.qTopP;
            if verify_conservation
                computeEnergyH(obj);
                obj.theta{kn_, 15} = obj.energy;
            end
            % save if last time step
            if k_ == length(obj.Fo)
                saveTheta(obj, k_);
            end
        end 
        function patchThetaK(obj, k_, thetaP)
            % compiles and saves data cell that contains results for the
            % simulation that matches Kevin's model 
            rW_ = []; thetaW_ = []; M = length(obj.rbarW);
            n = length(obj.zbarW); m = cell(M, 1); thetaWC_ = cell(M, 1);
            for i = 1:M, m{i} = length(obj.rbarW{i}); end
            for i = 1:length(obj.thetaW)
                rW_ = [rW_, obj.rbarW{i}];
%                 thetaWC_{i} = obj.thetaW{i}((k_-1)*n+1:k_*n, ...
%                                                  (k_-1)*m{i}+1:k_*m{i});
                thetaW_ = [thetaW_, obj.thetaW{i}];
%                 thetaW_ = [thetaW_, thetaWC_{i}];
            end
            n = length(obj.zbarW);
            mp = length(obj.rbarH); mw = length(rW_); m = mp + mw;
            theta_ = NaN*ones(n, m);
            % save and reinitialize theta if no more space
            if mod(k_-1, obj.ls) == 0 && k_ ~= 1
                saveThetaK(obj, k_-1);
                obj.theta = {};                
            end
            % normalize time step
            kn_ = mod(k_-1, obj.ls) + 1;
            % include particle/air temperature
            theta_(1:n, 1:mp) = thetaP;
            % include wall temperature
            theta_(1:n, mp+1:end) = thetaW_;
            % store data
            obj.thetaK{kn_, 1} = theta_;
            obj.thetaK{kn_, 2} = obj.zbarW;
            obj.thetaK{kn_, 3} = [obj.rbarH, rW_]; 
            obj.thetaK{kn_, 4} = k_;
            obj.thetaK{kn_, 5} = obj.Fo(k_);
            obj.thetaK{kn_, 6} = obj.thetaW;
            obj.thetaK{kn_, 7} = obj.qWall;
            obj.thetaK{kn_, 8} = obj.thetaA;
            obj.thetaK{kn_, 9} = obj.qLossW;
            % save if last time step
            if k_ == length(obj.Fo)
                saveThetaK(obj, k_);
            end
        end
        function thetaI_ = computeThetaI(obj, thetaT_, thetaC_, rtop_)
            % computes the average temperature of the top and center flow
            % channel boundaries to set the offset temperature that avoids
            % boundary osscilations in the Green's function
            tr = bulkTempR(obj, mean(thetaT_, 1), rtop_);
            tz = bulkTempZ(obj, mean(thetaC_, 2));
            thetaI_ = mean([tr, tz]);
        end
        function savePoints(obj)
            % saves temperature value at discrete points defined in zStore
            % and rStore
            if isempty(obj.thetaStore)                
                for i = 1:length(obj.zStore)
                    c = struct('z', obj.zStore(i), 'r', obj.rStore(i), ...
                               'Fo', obj.Fo, 'theta', []);
                    obj.thetaStore{i} = c;
                end
            end
            for i = 1:length(obj.zStore)
                obj.thetaStore{i}.theta = [obj.thetaStore{i}.theta, ...
                    obj.thetaZR(obj.zStore(i), obj.rStore(i))];
            end
        end
        function x = thetaZR(obj, z, r)
            % extracts current nondimensional temperature at (z, r) from
            % the thetaS, thetaT and thetaC arrays
            if r < obj.a0
                % pull from thetaC
                [~, i] = min(abs(obj.zcenter - z));
                [~, j] = min(abs(obj.rhat - r));
                x = obj.thetaC(end-i+1, j);
            elseif z > obj.ztop
                % pull from thetaT
                [~, i] = min(abs(obj.zhat + obj.ztop - z));
                [~, j] = min(abs(obj.rtop - r));
                x = obj.thetaT(i, j);
            elseif z > obj.ztop + obj.h
                % set equal to ambient temp inside storage bin
                x = obj.thetaA;
            else
                % pull from thetaS
                [~, i] = min(abs(obj.zbar - z));
                [~, j] = min(abs(obj.rbar - r));
                x = obj.thetaS(i, j);
            end
        end
        function assembleUz(obj)
            % sets z-velocity component for every cell in mesh
            if isempty(obj.wbar)
                computeWbar(obj);
            end
            obj.uz = cell(7, 1);
            % stagnant region
            zs = 0:dz_:obj.ztop; rs = obj.a0:dr_:obj.b;
            obj.uz{2} = zeros(length(zs), length(rs));
            % top flowing surface
            zt = obj.ztop:dz_:obj.ztop+obj.h; rt = obj.a:dr_:obj.b;
            obj.uz{3} = zeros(length(zt)-1, length(rt));
            % center flowing channel
            zc = 0:dz_:obj.ztop; rc = 1e-6:dr_:obj.a0;
            wbar_ = interp1(obj.zcenter, obj.wbar, zc)';
            obj.uz{4} = repmat(wbar_, 1, length(rc)-1);
            % mixing region
            zm = obj.ztop:dz_:obj.ztop+obj.h; rm = 1e-6:dr_:obj.a;
            obj.uz{5} = zeros(length(zm)-1, length(rm)-1);
            % full domain
            obj.uz{1} = [obj.uz{5}, obj.uz{3}; obj.uz{4}, obj.uz{2}];
            obj.uz{6} = [zs, zt(2:end)]; obj.uz{7} = [rc(1:end-1), rs];
        end
        function assembleUr(obj)
            % sets r-velocity component for every cell in mesh
            if isempty(obj.ubar)
                computeUbar(obj);
            end
            obj.ur = cell(7, 1);
            dz_ = obj.dz; dr_ = obj.dr;
            % stagnant region
            zs = 0:dz_:obj.ztop; rs = obj.a0:dr_:obj.b;
            obj.ur{2} = zeros(length(zs), length(rs));
            % top flowing region
            zt = obj.ztop:dz_:obj.ztop+obj.h; rt = obj.a:dr_:obj.b;
            ubar_ = interp1(obj.rtop, obj.ubar, rt);
            obj.ur{3} = repmat(ubar_, length(zt)-1, 1);
            % center flowing channel
            zc = 0:dz_:obj.ztop; rc = 1e-6:dr_:obj.a0;
            obj.ur{4} = zeros(length(zc), length(rc)-1);
            % mixing region
            zm = obj.ztop:dz_:obj.ztop+obj.h; rm = 1e-6:dr_:obj.a;
            obj.ur{5} = zeros(length(zm)-1, length(rm)-1);
            % full domain
            obj.ur{1} = [obj.ur{5}, obj.ur{3}; obj.ur{4}, obj.ur{2}]; 
        end
        function assembleEnergyD(obj)
            % sets non-dimensional temperature values and velocities 
            % in a matching mesh for computing the energy equation in 
            % discharge mode
            if isempty(obj.thetaMatchD), obj.thetaMatchD = cell(8, 1); end
            if isempty(obj.ubar), computeUbar(obj); end
            if isempty(obj.wbar), computeWbar(obj); end
            obj.ur = cell(5, 1);
            obj.uz = cell(5, 1);
            dz_ = obj.dz; dr_ = obj.dr;
            % stagnant region
            zs = 0:dz_:obj.ztop; rs = obj.a0:dr_:obj.b;
            [Rs, Zs] = meshgrid(obj.rbar, obj.zbar);
            [Rsq, Zsq] = meshgrid(rs, zs);
            obj.thetaMatchD{5} = interp2(Rs, Zs, obj.thetaS, ...
                                                   Rsq, Zsq, 'makima');
            obj.ur{2} = zeros(length(zs), length(rs));
            obj.uz{2} = zeros(length(zs), length(rs));
            % top flowing region
            zt = obj.ztop:dz_:obj.ztop+obj.h; rt = obj.a:dr_:obj.b;
            [Rt, Zt] = meshgrid(obj.rtop, obj.zhat);
            [Rtq, Ztq] = meshgrid(rt, zt(2:end));
            obj.thetaMatchD{6} = interp2(Rt, Zt, obj.thetaT, ...
                                                   Rtq, Ztq, 'makima');
            ubar_ = interp1(obj.rtop, obj.ubar, rt);
            obj.ur{3} = repmat(ubar_, length(zt)-1, 1);
            obj.uz{3} = zeros(length(zt)-1, length(rt));
            % center flowing channel
            zc = 0:dz_:obj.ztop; rc = 1e-6:dr_:obj.a0;
            [Rc, Zc] = meshgrid(obj.rhat, obj.zcenter);
            [Rcq, Zcq] = meshgrid(rc(1:end-1), zc);
            obj.thetaMatchD{7} = interp2(Rc, Zc, obj.thetaC, ...
                                                 Rcq, Zcq, 'makima');
            obj.ur{4} = zeros(length(zc), length(rc)-1);
            wbar_ = interp1(obj.zcenter, obj.wbar, zc)';
            obj.uz{4} = repmat(wbar_, 1, length(rc)-1);
            % mixing region
            zm = obj.ztop:dz_:obj.ztop+obj.h; rm = 1e-6:dr_:obj.a;
            [Rmq, Zmq] = meshgrid(rm(1:end-1), zm(2:end));
            obj.thetaMatchD{8} = obj.thetaChat*ones(length(zm(1:end-1)), ...
                                                        length(rm(2:end)));
            obj.ur{5} = zeros(length(zm)-1, length(rm)-1);
            obj.uz{5} = zeros(length(zm)-1, length(rm)-1);
            % full domain
            Z = obj.thetaMatchD{3}; R = obj.thetaMatchD{4};
            Zq = [Zcq, Zsq; Zmq, Ztq]; Rq = [Rcq, Rsq; Rmq, Rtq];
            if isempty(obj.thetaMatchD{2})                
                obj.thetaMatchD{1} = [obj.thetaMatchD{7}, obj.thetaMatchD{5}; ...
                                     obj.thetaMatchD{8}, obj.thetaMatchD{6}]; 
                obj.thetaMatchD{2} = obj.thetaMatchD{1};
            else
                obj.thetaMatchD{2} = interp2(R, Z, obj.thetaMatchD{1}, ...
                                                        Rq, Zq, 'makima');
                obj.thetaMatchD{1} = [obj.thetaMatchD{7}, obj.thetaMatchD{5}; ...
                                     obj.thetaMatchD{8}, obj.thetaMatchD{6}];             
            end 
            obj.ur{1} = [obj.ur{4}, obj.ur{2}; obj.ur{5}, obj.ur{3}];
            obj.uz{1} = [obj.uz{4}, obj.uz{2}; obj.uz{5}, obj.uz{3}];
            obj.thetaMatchD{3} = Zq; obj.thetaMatchD{4} = Rq;
        end
        function assembleEnergyH(obj)
            % sets non-dimensional temperature values in a matching mesh 
            % for computing the energy equation in charge and holding mode
            if isempty(obj.thetaMatchH), obj.thetaMatchH = cell(10, 1); end
            zp = 0:obj.dz:obj.ztop; rp = 1e-6:obj.dr:obj.b;
            [Rp, Zp] = meshgrid(obj.rbarH, obj.zbarH);
            [Rpq, Zpq] = meshgrid(rp, zp);
            t = interp2(Rp, Zp, full(obj.thetaH), Rpq, Zpq);
            if isempty(obj.thetaMatchH{2})
                obj.thetaMatchH{1} = t; obj.thetaMatchH{2} = t;
            else
                Z = obj.thetaMatchH{3}; R = obj.thetaMatchH{4}; 
                obj.thetaMatchH{2} = interp2(R, Z, obj.thetaMatchH{1}, ...
                                                      Rpq, Zpq, 'makima');
                obj.thetaMatchH{1} = t;
            end
            obj.thetaMatchH{3} = Zpq; obj.thetaMatchH{4} = Rpq;  
        end
        function [theta_, z_, r_] = assembleTheta(obj)
            % sets non-dimensional temperature values in a matching mesh
            % set matching mesh size as smallest dz, dr
%             dz_ = min([obj.dzbar, obj.dzc, obj.dzhat]);
            dz_ = min([obj.dzc, obj.dzhat]);
%             dr_ = min([obj.drbar, obj.drhat, obj.drtop]);
            dr_ = min([obj.drhat, obj.drtop]);
            % stagnant region
            zs = 0:dz_:obj.ztop; rs = obj.a0:dr_:obj.b;
            [Rs, Zs] = meshgrid(obj.rbar, obj.zbar);
            [Rsq, Zsq] = meshgrid(rs, zs);
            thetaS_ = interp2(Rs, Zs, obj.thetaS, Rsq, Zsq);
            % top flowing region
            zt = obj.ztop:dz_:obj.ztop+obj.h; rt = obj.a:dr_:obj.b;
            [Rt, Zt] = meshgrid(obj.rtop, obj.zhat);
            [Rtq, Ztq] = meshgrid(rt, zt(2:end));
            thetaT_ = interp2(Rt, Zt, obj.thetaT, Rtq, Ztq);
            % center flowing channel
            zc = 0:dz_:obj.ztop; rc = 1e-6:dr_:obj.a0;
            [Rc, Zc] = meshgrid(obj.rhat, obj.zcenter);
            [Rcq, Zcq] = meshgrid(rc(1:end-1), zc);
            thetaC_ = interp2(Rc, Zc, obj.thetaC, Rcq, Zcq);
            % mixing region (set to zero to null cell-to-cell analysis
            zm = obj.ztop:dz_:obj.ztop+obj.h; rm = 1e-6:dr_:obj.a;
            thetaM_ = obj.thetaChat*ones(length(zm(1:end-1)), ...
                                                        length(rm(2:end)));
            % full domain
            theta_ = [thetaC_, thetaS_; thetaM_, thetaT_];
            z_ = [zs, zt]; r_ = [rc, rs];
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % conversions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function T = theta2T(obj, theta)
            T = theta*(obj.T0 - obj.Tinf) + obj.Tinf;
        end
        function theta_ = T2theta(obj, T)
            theta_ = (T - obj.Tinf)/(obj.T0 - obj.Tinf);
        end
        function t = Fo2t(obj, Fo, prototype)
            if nargin < 3, prototype = 0; end
            if prototype 
                t = obj.Hp^2*Fo/obj.alphapPacked;
            else
                t = obj.H_^2*Fo/obj.alphaPacked;
            end
        end
        function Fo_ = t2Fo(obj, t, prototype)
            if nargin < 3, prototype = 0; end
            if prototype
                Fo_ = t*obj.alphapPacked/obj.Hp^2;
            else
                Fo_ = t*obj.alphaPacked/obj.H_^2;
            end
                
        end
        function y = S2T(~, thetaS_, r_, rtop_)
            % matches size of a row of stagnant region with the the size
            % of a row in the top boundary
            y = interp1(r_, thetaS_(end-1, :), rtop_, ...
                 'linear', 'extrap');
        end
        function y = S2C(~, thetaS_, z_, zcenter_)
            % matches size of a row of stagnant region with the the size
            % of a row in the top boundary
            y = interp1(z_, flipud(thetaS_(:, 1)), zcenter_, ...
                 'linear', 'extrap');
        end
        function y = T2S(~, thetaT_, rtop_, r_, averaging)
            % matches size of a row in the top boundary with the size of a
            % row in the stagnant region
            if nargin < 5
                averaging = 'full';
            end
            switch averaging
                case 'full'
                    y = NaN*ones(size(thetaT_, 1), length(r_));
                    for i = 1:size(y, 1)
                        y(i, :) = interp1(rtop_, thetaT_(i, :), ...
                            r_, 'linear', 'extrap');
                    end
                case 'mid'
                    y = interp1(rtop_, mean(thetaT_), ...
                        r_, 'linear', 'extrap');
                case 'top'
                    y = interp1(rtop_, thetaT_(1, :), ...
                        r_, 'linear', 'extrap');
                case 'bottom'
                    y = interp1(rtop_, thetaT_(end, :), ...
                        r_, 'linear', 'extrap');
            end
        end
        function y = C2S(~, thetaC_, zcenter_, z_, averaging)
            % matches size of a column in the center boundary with the 
            % size of a column in the stagnant region
            if nargin < 5
                averaging = 'full';
            end
            switch averaging
                case 'full'
                    y = NaN*ones(length(z_), size(thetaC_, 2));
                    for j = 1:size(y, 2)
                            y(:, j) = interp1(zcenter_, thetaC_(:, j), ...
                                z_, 'linear', 'extrap');
                    end
                case 'mid'
                    y = interp1(zcenter_, flipud(mean(thetaC_, 2)), ...
                        z_, 'linear', 'extrap');
                case 'top'
                    y = interp1(zcenter_, flipud(thetaC_(:, 1)), ...
                        z_, 'linear', 'extrap');
                case 'bottom'
                    y = interp1(zcenter_, flipud(thetaC_(:, end)), ...
                        z_, 'linear', 'extrap');
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % other
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
        function x = thetak(obj, k_)
            % returns theta matrix at time step k_
            ns = length(obj.zhat) - 1; nl = length(obj.zbar0); n = ns + nl;
            ms = length(obj.rhat) - 1; ml = length(obj.rbar); m = ms + ml;
            x = obj.theta(n*(k_-1)+1:n*k_, m*(k_-1)+1:m*k_);
        end       
        function saveTheta(obj, k_)
            % saves theta matrix at time step k_
            if k_ == length(obj.Fo)
                thetaString = sprintf('theta_%d.mat', floor(k_/obj.ls) + 1);
            else
                thetaString = sprintf('theta_%d.mat', floor(k_/obj.ls));
            end
            saveString = strcat(obj.thetaFolder, '\', thetaString);
            thetaSave = obj.theta;
            save(saveString, 'thetaSave');  
        end
        function saveThetaK(obj, k_)
            % saves theta matrix at time step k_
            if k_ == length(obj.Fo)
                thetaString = sprintf('thetaK_%d.mat', floor(k_/obj.ls) + 1);
            else
                thetaString = sprintf('thetaK_%d.mat', floor(k_/obj.ls));
            end
            saveString = strcat(obj.thetaFolder, '\', thetaString);
            thetaSave = obj.thetaK;
            save(saveString, 'thetaSave');  
        end
        function loadTheta(obj, k_)
            % loads theta matrix from stored location. The range of data
            % contained in the loaded matrix ends at time step k_
            if k_ == length(obj.Fo)
                thetaFile = sprintf('theta_%d.mat', floor(k_/obj.ls) + 1);
            else
                thetaFile = sprintf('theta_%d.mat', floor(k_/obj.ls));
            end
            thetaPath = strcat(obj.thetaFolder, '\', thetaFile);
            load(thetaPath, 'thetaSave');
            obj.theta = thetaSave;
        end 
        function loadThetaK(obj, k_)
            % loads theta matrix from stored location. The range of data
            % contained in the loaded matrix ends at time step k_
            if k_ == length(obj.Fo)
                thetaFile = sprintf('thetaK_%d.mat', floor(k_/obj.ls) + 1);
            else
                thetaFile = sprintf('thetaK_%d.mat', floor(k_/obj.ls));
            end
            thetaPath = strcat(obj.thetaFolder, '\', thetaFile);
            load(thetaPath, 'thetaSave');
            obj.thetaK = thetaSave;
        end 
        function tb = bulkTempR(obj, t, r_)
            % computes the bulk temperature for an r-dimensional array
            fr = simpsonIntegrator(obj, r_);
            tb = (t.*r_)*fr'/(r_*fr');
        end
        function tb = bulkTempZ(~, t)
            % computes the bulk temperature for an r-dimensional array
            tb = mean(t, 1);
        end
        function c = variableC(~, t)
            % computes the temperature-dependent specific heat for every
            % value in t. Outputs specific heat matrix of same dim as t. t
            % must be in Kelvin.
            a_ = 148.2; b_ = 0.3093; c = a_*t.^b_;
        end
        function w = simpsonWeights(~, p)
            % generates a set of weights from an array of points for the
            % 3-point quadrature rule used in the Simpson Integrator
            w = NaN*ones(size(p)); n = size(p, 1);
            % build weighting matrix neglecting first row
            a_ = p(:, 1); b_ = p(:, 2); c_ = p(:, 3);
            det = 1./(a_.^2.*(c_-b_) + b_.^2.*(a_-c_) + c_.^2.*(b_-a_));
            w(:, 1) = det.*((b_ - c_).*(b_.^3 - c_.^3)/3 + ...
                            (b_ - c_).*(b_.^2.*c_ - b_.*c_.^2) + ...
                            (c_.^2 - b_.^2).*(b_.^2 - c_.^2)/2);
            w(:, 2) = det.*((a_.^2 - c_.^2).*(b_.^2 - c_.^2)/2 + ...
                            (b_ - c_).*(a_.*c_.^2 - a_.^2.*c_) + ...
                            (c_ - a_).*(b_.^3 - c_.^3)/3);
            w(:, 3) = det.*((a_ - b_).*(b_.^3 - c_.^3)/3 + ...
                            (b_ - c_).*(a_.^2.*b_ - a_.*b_.^2) + ...
                            (b_.^2 - a_.^2).*(b_.^2 - c_.^2)/2);
            w = spdiags(w, [0, 1, 2], n, n+2);
            % add first row with reversed quadrature
            w = [zeros(1, n+2); w]; 
            w(1, 1) = w(2, 3); w(1, 2) = w(2, 2); w(1, 3) = w(2, 1); 
        end
        function x = simpsonIntegrator(obj, v)
            % generates a vector that can be used to integrate a matrix or
            % vector along the dimension with a spatial mesh given by v
            n = length(v); v = reshape(v, n, 1);
            p_ = [v(1:end-2), v(2:end-1), v(3:end)];
            w = simpsonWeights(obj, p_);  
            x = sum(w);
        end
        function [x, dv] = simpsonIntegratorEq(~, v)
            % generates a vector that can be used to integrate a matrix or
            % vector along the dimension with a spatial mesh given by v
            dv = abs(v(2) - v(1));
            x = NaN*ones(size(v));
            x(1) = dv/3; x(2) = 15*dv/12; x(3) = 11*dv/12;
            x(4:end-2) = dv; x(end-1) = 13*dv/12; x(end) = 5*dv/12;
        end
        function [x, dx] = nodeGen(~, xlim, n)
            % generates a set of n chebyshev nodes spaced between xlim(1)
            % and xlim(2)
            r_ = (xlim(2) - xlim(1))/2; theta_ = linspace(pi, 0, n);
            x = xlim(1) + r_*(1 + cos(theta_));
            dx = (eye(n) - diag(ones(1, n-1), -1))*x'; dx = dx(2:end);
        end
        function D = diffR(~, n, m, dr_)
            % creates a differential operator for first derivative in r 
            % dimension with central differencing at internal points 
            % and forward/backward differencing at boundaries
            nm = n*m;
            % set central differencing for center points
            D = spdiags([ones(nm, 1), -ones(nm, 1)], [1, -1], nm, nm);
            D = D./(2*dr_);
            % delete boundary entries
            D(m:m:end, :) = 0;      % r = 0
            D(m+1:m:end, :) = 0;    % r = b            
            % set forward differencing for left boundary
            D(1:m*(nm+1):end) = -1/dr_;
            D(nm+1:m*(nm+1):end) = 1/dr_;           
            % set backward differencing for right boundary
            D((m-1)*(nm+1)+1:m*(nm+1):end) = 1/dr_;
            D((m-2)*(nm+1)+2:m*(nm+1):end) = -1/dr_; 
        end   
        function D = diff2R(~, n, m, dr_)
            % creates a differential operator for first derivative in r 
            % dimension with central differencing at internal points 
            % and forward/backward differencing at boundaries
            nm = n*m;            
            % set central differencing for center points
            D = spdiags([-2*ones(nm, 1), ones(nm, 1), ones(nm, 1)], ...
                             [0, 1, -1], nm, nm);
            D = D./dr_^2;            
            % delete boundary entries
            D(m:m:end, :) = 0;      % r = 0
            D(m+1:m:end, :) = 0;    % r = b            
            % set forward differencing for left boundary
            D(1:m*(nm+1):end) = 1/dr_^2;
            D(nm+1:m*(nm+1):end) = -2/dr_^2;
            D(2*nm+1:m*(nm+1):end) = 1/dr_^2; 
            % set backward differencing for right boundary
            D((m-1)*(nm+1)+1:m*(nm+1):end) = 1/dr_^2;
            D((m-2)*(nm+1)+2:m*(nm+1):end) = -2/dr_^2;
            D((m-3)*(nm+1)+3:m*(nm+1):end) = 1/dr_^2; 
        end
        function D = diffZ(~, n, m, dz_)
            % creates a differential operator for first derivative in z 
            % dimension with central differencing at internal points 
            % and forward/backward differencing at boundaries
            nm = n*m;            
            % set central differencing for center points
            D = spdiags([ones(nm, 1), -ones(nm, 1)], [m, -m], nm, nm);
            D = D./(2*dz_);            
            % delete boundary entries
            D(1:m, :) = 0;          % z = 0
            D(end-m+1:end, :) = 0;    % z = 1            
            % set forward differencing for top boundary
            D(1:nm+1:m*(nm+1)) = -1/dz_;
            D(m*nm+1:nm+1:2*m*nm+1) = 1/dz_;            
            % set backward differencing for bottom boundary
            D(end-(m-1)*(nm+1):nm+1:end) = 1/dz_;
            D(end-2*m*(nm+1)+nm+1+m:nm+1:end-m*nm) = -1/dz_;               
        end
        function D = diff2Z(~, n, m, dz_)
            % creates a differential operator for first derivative in z 
            % dimension with central differencing at internal points 
            % and forward/backward differencing at boundaries
            nm = n*m;            
            % set central differencing for center points
            D = spdiags([-2*ones(nm, 1), ones(nm, 1), ones(nm, 1)], ...
                          [0, m, -m], nm, nm);
            D = D./dz_^2;            
            % delete boundary entries
            D(1:m, :) = 0;            % z = 0
            D(end-m+1:end, :) = 0;    % z = 1            
            % set forward differencing for top boundary
            D(1:nm+1:m*(nm+1)) = 1/dz_^2;
            D(m*nm+1:nm+1:2*m*nm+1) = -2/dz_^2;
            D(2*m*nm+1:nm+1:3*m*nm+1) = 1/dz_^2;   
            % set backward differencing for bottom boundary
            D(end-(m-1)*(nm+1):nm+1:end) = 1/dz_^2;
            D(end-2*m*(nm+1)+nm+1+m:nm+1:end-m*nm) = -2/dz_^2;
            D(end-3*m*(nm+1)+nm+1+2*m:nm+1:end-2*m*nm) = 1/dz_^2;
        end      
        % static analytic solutions         
        function t = Xt(~, Fo_, beta_, eta_)
            % transient component of analytic solutions
            t = exp(-(beta_^2 + eta_^2)*Fo_);
        end
    end
end

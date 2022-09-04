model FallingParticleReceiverSystem
  package Types
    extends Icons.Library;
    type MassFlowRate = Real(final quantity = "MassFlowRate", final unit = "kg/s");
    type Temperature = Real(final quantity = "ThermodynamicTemperature", final unit = "K", min = 0.0, start = 288.15, nominal = 300, displayUnit = "degC") "Absolute temperature (use type TemperatureDifference for relative temperatures)" annotation(
      absoluteValue = true);
    type SpecificHeatCapacity = Real(final quantity = "SpecificHeatCapacity", final unit = "J/(kg.K)");
    type Velocity = Real(final quantity = "Velocity", final unit = "m/s");
    type Length = Real(final quantity = "Length", final unit = "m");
    type Heat = Real(final quantity = "Heat", final unit = "W");
    type HeatFlux = Real(final quantity = "HeatFlux", final unit = "W/m2");
    type Area = Real(final quantity = "Area", final unit = "m2");
    type Enthalpy = Real(final quantity = "Enthalpy", final unit = "J/kg");
    type Volume = Real(final quantity = "Volume", final unit = "m3");
    type Density = Real(final quantity = "Density", final unit = "kg/m3");
    type Convection = Real(final quantity = "ConvectionCoefficient", final unit = "W/m2-K");
    type Mass = Real(final quantity = "Mass", final unit = "kg");
    type Fraction = Real(final quantity = "Fraction", final unit = "N/A");
    type Energy = Real(final quantity = "Energy", final unit = "MW.hr");
    type Conductivity = Real(final quantity = "ThermalConductivity", final unit = "W/m-K");
    type Day = Real(final quantity = "Day", final unit = "day");
    type Hour = Real(final quantity = "Hour", final unit = "hr");
    type Second = Real(final quantity = "Second", final unit = "s");
  end Types;

  package Interfaces
    extends Icons.Library;

    connector ParticleFlow
      flow Types.MassFlowRate m_dot;
      Types.Temperature T;
      annotation(
        Icon(graphics = {Ellipse(fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}}, endAngle = 360)}));
    end ParticleFlow;

    connector Temperature
      Types.Temperature T;
      annotation(
        Icon(graphics = {Ellipse(fillColor = {255, 85, 0}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}}, endAngle = 360)}));
    end Temperature;

    connector CO2Flow
      flow Types.MassFlowRate m_dot;
      Types.Temperature T;
      annotation(
        Icon(graphics = {Ellipse(fillColor = {0, 0, 255}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}}, endAngle = 360)}));
    end CO2Flow;

    connector MassFlow
      Types.MassFlowRate m_dot annotation(
        Diagram,
        Icon(graphics = {Ellipse(fillColor = {0, 85, 0}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}}, endAngle = 360)}));
      annotation(
        Icon(graphics = {Ellipse(fillColor = {0, 85, 0}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}}, endAngle = 360)}));
    end MassFlow;

    connector Heat
      Types.Heat Q;
      annotation(
        Icon(graphics = {Ellipse(fillColor = {255, 0, 0}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}}, endAngle = 360)}, coordinateSystem(initialScale = 0.1)));
    end Heat;
  end Interfaces;

  package SourceSink
    extends Icons.Library;

    model ParticleSource
      extends Icons.SolidRes;
      Types.Temperature T;
      Types.MassFlowRate m_dot;
      Interfaces.ParticleFlow ParticleOutlet annotation(
        Placement(visible = true, transformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      T = ParticleOutlet.T;
      m_dot = -ParticleOutlet.m_dot;
    end ParticleSource;

    model ParticleSink
      extends Icons.SolidRes;
      Types.Temperature T;
      Types.MassFlowRate m_dot;
      FallingParticleReceiverSystem.Interfaces.ParticleFlow ParticleInlet annotation(
        Placement(visible = true, transformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      T = ParticleInlet.T;
      m_dot = ParticleInlet.m_dot;
    end ParticleSink;

    model CO2Source
      Types.Temperature T;
      Types.MassFlowRate m_dot;
      FallingParticleReceiverSystem.Interfaces.CO2Flow CO2Outlet annotation(
        Placement(visible = true, transformation(origin = {100, -2}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      m_dot = -CO2Outlet.m_dot;
      T = CO2Outlet.T;
      annotation(
        Icon(graphics = {Polygon(fillColor = {0, 0, 255}, fillPattern = FillPattern.Solid, points = {{-100, 60}, {-100, -60}, {40, -60}, {100, 0}, {40, 60}, {-100, 60}})}));
    end CO2Source;

    model CO2Sink
      Types.Temperature T;
      Types.MassFlowRate m_dot;
      FallingParticleReceiverSystem.Interfaces.CO2Flow CO2Inlet annotation(
        Placement(visible = true, transformation(origin = {-98, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      T = CO2Inlet.T;
      m_dot = CO2Inlet.m_dot;
      annotation(
        Icon(graphics = {Polygon(fillColor = {0, 0, 255}, fillPattern = FillPattern.Solid, points = {{-100, 0}, {-40, -60}, {100, -60}, {100, 60}, {-40, 60}, {-40, 60}, {-100, 0}})}));
    end CO2Sink;
  end SourceSink;

  package Media
    extends Icons.Library;

    package Particle
      extends Icons.PropertyPackage;

      function Enthalpy
        extends Icons.Function;
        input Types.Temperature T;
        output Types.Enthalpy h;
      protected
        parameter Types.SpecificHeatCapacity cp_s = 1250;
        parameter Real T_ref = 273.15;
      algorithm
        h := cp_s * (T - T_ref);
      end Enthalpy;

      function Density
        extends Icons.Function;
        output Types.Density rho;
      protected
        parameter Types.Density rho_s = 3500;
      algorithm
        rho := rho_s;
      end Density;

      function SpecificHeat
        extends Icons.Function;
        input Types.Temperature T;
        output Types.Enthalpy cp;
      protected
        parameter Types.SpecificHeatCapacity cp_s = 1250;
      algorithm
        cp := cp_s;
      end SpecificHeat;

      model PackingLimit
        extends Icons.Function;
        output Types.Density phi_s;
      protected
        parameter Types.Density rho_s = 0.60;
      algorithm
        phi_s := phi;
      end PackingLimit;
    end Particle;

    package CO2
      extends Icons.PropertyPackage;

      function Density
        extends Icons.Function;
        output Types.SpecificHeatCapacity rho;
      protected
        parameter Types.Density rho_CO2 = 110;
      algorithm
        rho := rho_CO2;
      end Density;

      function Viscosity
        extends Icons.Function;
      end Viscosity;

      function Conductivity
        extends Icons.Function;
        output Types.SpecificHeatCapacity k;
      protected
        parameter Types.SpecificHeatCapacity k_CO2 = 0.05;
      algorithm
        k := k_CO2;
      end Conductivity;

      function SpecificHeat
        extends Icons.Function;
        output Types.SpecificHeatCapacity cp;
      protected
        parameter Types.SpecificHeatCapacity cp_CO2 = 1250;
      algorithm
        cp := cp_CO2;
      end SpecificHeat;

      function Enthalpy
        extends Icons.Function;
        input Types.Temperature T;
        output Types.Enthalpy h;
      protected
        parameter Types.SpecificHeatCapacity cp_CO2 = 1250;
        parameter Real T_ref = 273.15;
      algorithm
        h := cp_CO2 * (T - T_ref);
      end Enthalpy;
    end CO2;

    package StainlessSteel
      extends Icons.PropertyPackage;

      function Density
        extends Icons.Function;
        output Types.Density rho;
      protected
        parameter Types.Density rho_m = 8000;
      algorithm
        rho := rho_m;
      end Density;

      function Conductivity
        extends Icons.Function;
        output Types.Conductivity k;
      protected
        parameter Types.Conductivity k_m = 0.5;
      algorithm
        k := k_m;
      end Conductivity;

      function SpecificHeat
        extends Icons.Function;
        output Types.SpecificHeatCapacity cp;
      protected
        parameter Types.SpecificHeatCapacity cp_m = 500;
      algorithm
        cp := cp_m;
      end SpecificHeat;
    end StainlessSteel;

    package RSLE_57
      extends Icons.PropertyPackage;

      function Density
        extends Icons.Function;
        output Types.Density rho;
      protected
        parameter Types.Density rho_m = 2000;
      algorithm
        rho := rho_m;
      end Density;

      function Conductivity
        extends Icons.Function;
        output Types.Conductivity k;
      protected
        parameter Types.Conductivity k_m = 0.5;
      algorithm
        k := k_m;
      end Conductivity;

      function SpecificHeat
        extends Icons.Function;
        output Types.SpecificHeatCapacity cp;
      protected
        parameter Types.SpecificHeatCapacity cp_m = 500;
      algorithm
        cp := cp_m;
      end SpecificHeat;
    end RSLE_57;

    package Superwool
      extends Icons.PropertyPackage;

      model Density
        extends Icons.Function;
      end Density;

      model Conductivity
        extends Icons.Function;
      end Conductivity;

      model SpecificHeat
        extends Icons.Function;
      end SpecificHeat;
    end Superwool;
  end Media;

  package Components
    extends Icons.Library;

    model FallingParticleReceiver
      extends Icons.CavityReceiver;
      //Variables
      parameter Types.Temperature T_o = 20 + 273.15;
      parameter Types.Temperature T_s_0 = 20 + 273.15 "initial temperature";
      Types.MassFlowRate m_dot_s_in "inlet mass flow rate";
      Types.MassFlowRate m_dot_s_out "outlet mass flow rate";
      Types.Temperature T_s_in "inlet particle temperature";
      Types.Temperature T_s_out(start = T_s_0) "outlet particle temperature";
      Types.Enthalpy h_s_in "inlet particle enthalpy";
      Types.Enthalpy h_s_out "outlet particle enthalpy";
      Types.Density rho_s "density";
      Types.SpecificHeatCapacity cp_s "specific heat capacity";
      Types.Heat Q_solar "solar irradiance";
      Types.Heat Q_rerad "reradiation heat loss";
      Types.Heat Q_in "receiver net input";
      parameter Types.Area A_ap = 1 "aperture area";
      Types.Fraction eta_rec "receiver efficiency";
      parameter Types.Volume V_rec = 1 "receiver volume";
      parameter Types.Fraction sigma_sb = 5.67 * 10 ^ (-8) "stephan boltzman constant";
      parameter Types.Fraction eps_s = 0.88 "emissivity";
      parameter Types.Fraction abs_s = 0.92 "emissivity";
      FallingParticleReceiverSystem.Interfaces.ParticleFlow ParticleInlet annotation(
        Placement(visible = true, transformation(origin = {0, 80}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, 80}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Interfaces.ParticleFlow ParticleOutlet annotation(
        Placement(visible = true, transformation(origin = {0, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealOutput OutletTemperature annotation(
        Placement(visible = true, transformation(origin = {-80, -78}, extent = {{-10, -10}, {10, 10}}, rotation = 180), iconTransformation(origin = {-80, -60}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealInput MassFlow annotation(
        Placement(visible = true, transformation(origin = {-80, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-80, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
//Connections
      m_dot_s_in = ParticleInlet.m_dot;
      m_dot_s_out = -ParticleOutlet.m_dot;
      T_s_in = ParticleInlet.T;
      T_s_out = ParticleOutlet.T;
      T_s_out = OutletTemperature;
      m_dot_s_in = MassFlow;
//Mass balance
      m_dot_s_in = m_dot_s_out;
//Heat balance
      Q_rerad = A_ap * sigma_sb * eps_s * (T_s_out ^ 4 - T_o ^ 4);
      Q_in = Q_solar * abs_s - Q_rerad;
      V_rec * rho_s * cp_s * der(T_s_out) = m_dot_s_in * h_s_in + Q_solar * abs_s - m_dot_s_out * h_s_out - Q_rerad;
      eta_rec * (Q_solar * abs_s) = m_dot_s_out * h_s_out - m_dot_s_in * h_s_in;
//Properties
      h_s_in = Media.Particle.Enthalpy(T_s_in);
      h_s_out = Media.Particle.Enthalpy(T_s_out);
      cp_s = Media.Particle.SpecificHeat(T_s_out);
      rho_s = Media.Particle.Density();
    end FallingParticleReceiver;

    model SolidSplitter
      extends Icons.Tee;
      //Variables
      Types.MassFlowRate m_dot_s_in "inlet mass flow rate";
      Types.Temperature T_s_in "inlet particle temperature";
      Types.MassFlowRate m_dot_s_out_1 "stream one outlet mass flow rate";
      Types.MassFlowRate m_dot_s_out_2 "stream two outlet mass flow rate";
      Types.Temperature T_s_out_1 "stream one outlet temperature";
      Types.Temperature T_s_out_2 "stream two outlet temperature";
      FallingParticleReceiverSystem.Interfaces.ParticleFlow ParticleInlet annotation(
        Placement(visible = true, transformation(origin = {-78, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-98, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      FallingParticleReceiverSystem.Interfaces.ParticleFlow ParticleOutlet1 annotation(
        Placement(visible = true, transformation(origin = {80, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      FallingParticleReceiverSystem.Interfaces.ParticleFlow ParticleOutlet2 annotation(
        Placement(visible = true, transformation(origin = {80, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, -100}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
//Connections
      m_dot_s_out_1 = -ParticleOutlet1.m_dot;
      m_dot_s_out_2 = -ParticleOutlet2.m_dot;
      m_dot_s_in = ParticleInlet.m_dot;
      T_s_out_1 = ParticleOutlet1.T;
      T_s_out_2 = ParticleOutlet2.T;
      T_s_in = ParticleInlet.T;
//Mass Balance
      m_dot_s_in = m_dot_s_out_1 + m_dot_s_out_2;
//Energy Balance
      T_s_in = T_s_out_1;
      T_s_in = T_s_out_2;
    end SolidSplitter;

    model SolidJunction
      extends Icons.Tee;
      //Variables
      Types.MassFlowRate m_dot_s_in_1 "stream one inlet mass flow rate";
      Types.MassFlowRate m_dot_s_in_2 "stream two inlet mass flow rate";
      Types.Temperature T_s_in_1 "stream one inlet temperature";
      Types.Temperature T_s_in_2 "stream two inlet temperature";
      Types.Enthalpy h_s_in_1 "stream one inlet enthalpy";
      Types.Enthalpy h_s_in_2 "stream two inlet enthalpy";
      Types.MassFlowRate m_dot_s_out "outlet mass flow rate";
      Types.Temperature T_s_out "outlet temperature";
      Types.Enthalpy h_s_out "outlet enthalpy";
      FallingParticleReceiverSystem.Interfaces.ParticleFlow ParticleInlet1 annotation(
        Placement(visible = true, transformation(origin = {-80, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      FallingParticleReceiverSystem.Interfaces.ParticleFlow ParticleInlet2 annotation(
        Placement(visible = true, transformation(origin = {-80, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, -100}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      FallingParticleReceiverSystem.Interfaces.ParticleFlow ParticleOutlet annotation(
        Placement(visible = true, transformation(origin = {78, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
//Connections
      m_dot_s_in_1 = ParticleInlet1.m_dot;
      m_dot_s_in_2 = ParticleInlet2.m_dot;
      m_dot_s_out = -ParticleOutlet.m_dot;
      T_s_in_1 = ParticleInlet1.T;
      T_s_in_2 = ParticleInlet2.T;
      T_s_out = ParticleOutlet.T;
//Mass Balance
      m_dot_s_in_1 + m_dot_s_in_2 = m_dot_s_out;
//Energy Balance
      m_dot_s_in_1 * h_s_in_1 + m_dot_s_in_2 * h_s_in_2 = m_dot_s_out * h_s_out;
//Properties
      h_s_in_1 = Media.Particle.Enthalpy(T_s_in_1);
      h_s_in_2 = Media.Particle.Enthalpy(T_s_in_2);
      h_s_out = Media.Particle.Enthalpy(T_s_out);
    end SolidJunction;

    model LumpedStorageBin
      extends Icons.StorageBin;
      //Geometric parameters
      parameter Types.Temperature T_o = 20 + 273.15;
      parameter Types.Mass m_s_0 = 1 "initial mass inventory";
      parameter Types.Temperature T_s_0 = 1 "initial temperature";
      parameter Types.Volume V_bin = 1 "storage bin volume";
      parameter Types.Area A_bin = 1 "storage bin cross section";
      parameter Types.Length H_bin = V_bin / A_bin "storage bin height";
      parameter Types.Convection h_loss = 10;
      parameter Types.Area A_s = 2 * A_bin + H_bin * (2 * sqrt(A_bin / 3.14) * 3.14) "storage bin surface area";
      //Particle variables
      Types.MassFlowRate m_dot_s_in "inlet mass flow rate";
      Types.MassFlowRate m_dot_s_out "outlet mass flow rate";
      Types.Temperature T_s_in "inlet temperature";
      Types.Temperature T_s_out(start = T_s_0) "outlet temperature";
      Types.Enthalpy h_s_in "inlet enthalpy";
      Types.Enthalpy h_s_out "outlet enthalpy";
      Types.Density rho_s "density";
      Types.SpecificHeatCapacity cp_s "specific heat capacity";
      parameter Real phi_s = 0.6 "solids volume fraction";
      Types.Mass m_s(start = m_s_0) "storage bin mass inventory";
      Types.Fraction X_fill "fill fraction";
      Types.Heat Q_loss "heat loss";
      FallingParticleReceiverSystem.Interfaces.ParticleFlow ParticleOutlet annotation(
        Placement(visible = true, transformation(origin = {0, -100}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, -100}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      FallingParticleReceiverSystem.Interfaces.ParticleFlow ParticleInlet annotation(
        Placement(visible = true, transformation(origin = {0, 98}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, 100}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
//Connection
      m_dot_s_in = ParticleInlet.m_dot;
      m_dot_s_out = -ParticleOutlet.m_dot;
      T_s_out = ParticleOutlet.T;
      T_s_in = ParticleInlet.T;
//Mass Balance
      der(m_s) = m_dot_s_in - m_dot_s_out;
      X_fill = m_s / (rho_s * phi_s) / V_bin;
//Energy Balance
      h_s_out * der(m_s) + cp_s * m_s * der(T_s_out) = m_dot_s_in * h_s_in - m_dot_s_out * h_s_out - h_loss * A_s * (T_s_out - T_o);
      Q_loss = h_loss * A_s * (T_s_out - T_o);
//Properties
      h_s_in = Media.Particle.Enthalpy(T_s_in);
      h_s_out = Media.Particle.Enthalpy(T_s_out);
      cp_s = Media.Particle.SpecificHeat(T_s_out);
      rho_s = Media.Particle.Density();
    end LumpedStorageBin;

    model ParticleHeatExchanger
      extends Icons.HeatExchanger;
      //Geometry variables
      parameter Integer n = 100 "number of discretizations";
      parameter Types.Length hc_s = 0.006 "solid channel width";
      parameter Types.Length hc_CO2 = 0.001 "CO2 channel width";
      parameter Types.Length t_m = 0.003 "metal thickness";
      parameter Types.Length H = 1 "heat exchanger height";
      parameter Types.Length W = 0.5 "heat exchanger width";
      parameter Integer N_plate = 25 "number of parallel plates";
      parameter Types.Convection h_conv_CO2 = 1000 "CO2 convection coefficient";
      parameter Types.Convection h_conv_sw = 180 "solid-wall convection coefficient";
      parameter Real dx = H / (n - 2) "discretization size";
      //Particle variables
      Types.MassFlowRate m_dot_s_in "inlet mass flow rate";
      Types.MassFlowRate m_dot_s_out "outlet mass flow rate";
      Types.Temperature T_s_in "inlet temperature";
      Types.Temperature T_s_out "outlet temperature";
      Types.Temperature T_s[n](each start = 273.15 + 25) "solid temperature distribution";
      Types.Enthalpy h_s[n] "particle enthalpy";
      Types.SpecificHeatCapacity cp_s "particle specific heat";
      Types.Density rho_s "particle density";
      parameter Real phi_s = 0.6 "solid volume fraction";
      Types.Velocity v_s "solid velocity";
      //CO2 variables
      Types.MassFlowRate m_dot_CO2_in "inlet mass flow rate";
      Types.MassFlowRate m_dot_CO2_out "outlet mass flow rate";
      Types.Temperature T_CO2_in "inlet temperature";
      Types.Temperature T_CO2_out "outlet temperature";
      Types.Temperature T_CO2[n](each start = 273.15 + 25) "CO2 temperature distribution";
      Types.Enthalpy h_CO2[n] "CO2 enthalpy";
      Types.SpecificHeatCapacity cp_CO2 "CO2 specific heat";
      Types.Density rho_CO2 "CO2 density";
      Types.Velocity v_CO2 "CO2 velocity";
      Types.Heat Q_CO2 "sCO2 heat addition";
      //Heat transfer surface
      Types.Temperature T_m[n](each start = 273.15 + 25) "metal temperature distribution";
      Types.SpecificHeatCapacity cp_m "metal specific heat";
      Types.Density rho_m "metal density";
      //Connections
      FallingParticleReceiverSystem.Interfaces.ParticleFlow ParticleInlet annotation(
        Placement(visible = true, transformation(origin = {0, 100}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, 132}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      FallingParticleReceiverSystem.Interfaces.ParticleFlow ParticleOutlet annotation(
        Placement(visible = true, transformation(origin = {0, -100}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, -126}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      FallingParticleReceiverSystem.Interfaces.CO2Flow CO2Inlet annotation(
        Placement(visible = true, transformation(origin = {-98, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-138, 6}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      FallingParticleReceiverSystem.Interfaces.CO2Flow CO2Outlet annotation(
        Placement(visible = true, transformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {132, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealOutput CO2OutletTemperature annotation(
        Placement(visible = true, transformation(origin = {98, -22}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {136, -24}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealOutput ParticleOutletTemperature annotation(
        Placement(visible = true, transformation(origin = {-48, -88}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-28, -128}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      Modelica.Blocks.Interfaces.RealInput m_dot_CO2 annotation(
        Placement(visible = true, transformation(origin = {36, -88}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-26, 130}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    equation
//Connections
      m_dot_s_in = ParticleInlet.m_dot;
      m_dot_s_out = -ParticleOutlet.m_dot;
      T_s_in = ParticleInlet.T;
      T_s_out = ParticleOutlet.T;
      m_dot_CO2_in = CO2Inlet.m_dot;
      m_dot_CO2_in = m_dot_CO2;
      m_dot_CO2_out = -CO2Outlet.m_dot;
      T_CO2_in = CO2Inlet.T;
      T_CO2_out = CO2Outlet.T;
      T_CO2_out = CO2OutletTemperature;
      T_s_out = ParticleOutletTemperature;
//Mass Balance
      m_dot_s_in = m_dot_s_out; 
      m_dot_s_in = rho_s * phi_s * v_s * hc_s * N_plate * W;
      m_dot_CO2_in = m_dot_CO2_out;
      m_dot_CO2_in = rho_CO2 * v_CO2 * hc_CO2 * N_plate * W;
//Energy Balance
      T_s[1] = T_s_in;
      for i in 2:n - 1 loop
        cp_s * rho_s * phi_s * der(T_s[i]) = (-cp_s * rho_s * phi_s * v_s * (T_s[i] - T_s[i - 1]) / dx) + 2 * h_conv_sw / hc_s * (T_m[i] - T_s[i]);
      end for;
      T_s[n] = T_s[n - 1];
      T_s[n] = T_s_out;
      T_CO2[1] = T_CO2_out;
      T_CO2[1] = T_CO2[2];
      for i in 2:n - 1 loop
        cp_CO2 * rho_CO2 * der(T_CO2[i]) = (-cp_CO2 * rho_CO2 * v_CO2 * (T_CO2[i] - T_CO2[i + 1]) / dx) + 2 * h_conv_CO2 / hc_CO2 * (T_m[i] - T_CO2[i]);
      end for;
      T_CO2[n] = T_CO2_in;
      T_m[1] = T_m[2];
      for i in 2:n - 1 loop
        cp_m * rho_m * der(T_m[i]) = 2 * h_conv_CO2 / t_m * (T_CO2[i] - T_m[i]) + 2 * h_conv_sw / t_m * (T_s[i] - T_m[i]);
      end for;
      T_m[n] = T_m[n - 1];
      Q_CO2 = cp_CO2 * m_dot_CO2_in * (T_CO2_out - T_CO2_in);
//Properties
      for i in 1:n loop
        h_s[i] = Media.Particle.Enthalpy(T_s[i]);
        h_CO2[i] = Media.CO2.Enthalpy(T_CO2[i]);
      end for;
      cp_s = Media.Particle.SpecificHeat(T_s_out);
      rho_s = Media.Particle.Density();
      cp_CO2 = Media.CO2.SpecificHeat();
      rho_CO2 = Media.CO2.Density();
      cp_m = Media.StainlessSteel.SpecificHeat();
      rho_m = Media.StainlessSteel.Density();
    end ParticleHeatExchanger;

    model FreeFallDuct
      extends Icons.InsulatedDuct;
      parameter Integer n = 10 "number of discretizations";
      parameter Types.Temperature T_o = 20 + 273.15 "ambient temperature";
      parameter Types.Temperature T_m_0 = 20 + 273.15 "initial temperature";
      parameter Types.Length perimeter = 1.0 "duct perimeter";
      parameter Types.Length length = 1.0 "duct length";
      parameter Real dx = length / (n - 2) "discretization size";
      parameter Types.Length thickness = 1.0 "duct thickness";
      parameter Types.Area A_s = perimeter * length "duct surface area";
      parameter Types.Volume V_m = thickness * A_s "duct metal volume";
      parameter Types.Convection h_loss = 10.0 "external heat loss convetion coefficient";
      parameter Types.Convection h_s_w = 100 "solid-wall convection coefficient";
      Types.MassFlowRate m_dot_s_in "inlet mass flow rate";
      Types.MassFlowRate m_dot_s_out "outlet mass flow rate";
      Types.Temperature T_s_in "inlet temperature";
      Types.Temperature T_s_out "outlet temperature";
      Types.Temperature T_m(start = T_m_0) "duct temperature";
      Types.Enthalpy h_s_in "inlet enthalpy";
      Types.Enthalpy h_s_out "outlet enthalpy";
      Types.SpecificHeatCapacity cp_m "duct specific heat";
      Types.Density rho_m "duct density";
      Types.Heat Q_loss "heat loss";
      Interfaces.ParticleFlow ParticleInlet annotation(
        Placement(visible = true, transformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Interfaces.ParticleFlow ParticleOutlet annotation(
        Placement(visible = true, transformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
// Connections
      m_dot_s_in = ParticleInlet.m_dot;
      m_dot_s_out = -ParticleOutlet.m_dot;
      T_s_in = ParticleInlet.T;
      T_s_out = ParticleOutlet.T;
// Mass Balance
      m_dot_s_in = m_dot_s_out;
// Energy Balance
      cp_m * V_m * rho_m * der(T_m) = h_s_w * A_s * (T_s_out - T_m) + h_loss * A_s * (T_o - T_m);
      h_s_w * A_s * (T_s_out - T_m) + m_dot_s_out * h_s_out = m_dot_s_in * h_s_in;
      Q_loss = h_loss * A_s * (T_m - T_o);
// Properties
      h_s_in = Media.Particle.Enthalpy(T_s_in);
      h_s_out = Media.Particle.Enthalpy(T_s_out);
      cp_m = Media.StainlessSteel.SpecificHeat();
      rho_m = Media.StainlessSteel.Density();
    end FreeFallDuct;

    model WedgeMassFlowHopper
      extends Icons.MassFlowHopper;
      //Geometric Parameters
      parameter Integer n = 100;
      parameter Types.Volume V_bin = 1;
      parameter Types.Area A_bin = 1;
      parameter Types.Length H_bin = V_bin / A_bin;
      parameter Real dx = H_bin / (n - 2);
      parameter Types.Temperature T_s_0 = 20 + 273.15;
      //Particle Variables
      Types.MassFlowRate m_dot_s_in;
      Types.MassFlowRate m_dot_s_out;
      Types.Temperature T_s_in;
      Types.Temperature T_s[n](each start = T_s_0);
      Types.Temperature T_s_out;
      Types.Enthalpy h_s_in;
      Types.Enthalpy h_s_out;
      Types.Velocity v_s;
      Types.Density rho_s;
      Types.SpecificHeatCapacity cp_s;
      parameter Real phi_s = 0.6;
      FallingParticleReceiverSystem.Interfaces.ParticleFlow ParticleOutlet annotation(
        Placement(visible = true, transformation(origin = {0, -100}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, -100}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      FallingParticleReceiverSystem.Interfaces.ParticleFlow ParticleInlet annotation(
        Placement(visible = true, transformation(origin = {0, 98}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, 100}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
//Connection
      m_dot_s_in = ParticleInlet.m_dot;
      m_dot_s_out = -ParticleOutlet.m_dot;
      T_s_out = ParticleOutlet.T;
      T_s_in = ParticleInlet.T;
//Mass Balance
      m_dot_s_in = m_dot_s_out;
      m_dot_s_in = rho_s * phi_s * v_s * A_bin;
//Energy Balance
      T_s[1] = T_s_in;
      for i in 2:n loop
        cp_s * rho_s * phi_s * der(T_s[i]) = -v_s * rho_s * phi_s * cp_s * (T_s[i] - T_s[i - 1]) / dx;
      end for;
      T_s[n] = T_s_out;
//Property Evaluation
      h_s_in = Media.Particle.Enthalpy(T_s_in);
      h_s_out = Media.Particle.Enthalpy(T_s_out);
      cp_s = Media.Particle.SpecificHeat(T_s_out);
      rho_s = Media.Particle.Density();
    end WedgeMassFlowHopper;

    model ElectricalHeater
      extends Icons.Heater;
      //Particle variables
      Types.MassFlowRate m_dot_s_in "inlet mass flow rate";
      Types.MassFlowRate m_dot_s_out "outlet mass flow rate";
      Types.Temperature T_s_in "inlet temperature";
      Types.Temperature T_s_out "outlet temperature";
      Types.Enthalpy h_s_in "inlet enthalpy";
      Types.Enthalpy h_s_out "outlet enthalpy";
      Types.Heat Q_in "electrical heat addition";
      FallingParticleReceiverSystem.Interfaces.ParticleFlow ParticleOutlet annotation(
        Placement(visible = true, transformation(origin = {0, -100}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, -100}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      FallingParticleReceiverSystem.Interfaces.ParticleFlow ParticleInlet annotation(
        Placement(visible = true, transformation(origin = {0, 98}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, 100}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
//Connection
      m_dot_s_in = ParticleInlet.m_dot;
      m_dot_s_out = -ParticleOutlet.m_dot;
      T_s_out = ParticleOutlet.T;
      T_s_in = ParticleInlet.T;
//Mass Balance
      m_dot_s_in = m_dot_s_out;
//Energy Balance
      Q_in = m_dot_s_in * (h_s_out - h_s_in);
//Properties
      h_s_in = Media.Particle.Enthalpy(T_s_in);
      h_s_out = Media.Particle.Enthalpy(T_s_out);
    end ElectricalHeater;

    model BucketElevator
      extends Icons.BucketElevator;
      // Geometry
      parameter Integer n = 100 "number of discretizations";
      parameter Types.Temperature T_o = 20 + 273.15 "ambient temperature";
      parameter Types.Length w = 1.0 "elevator case width";
      parameter Types.Length l = 2.0 "elevator case length";
      parameter Types.Length h = 50.0 "elevator case height";
      parameter Types.Length t_m = 0.006 "elevator case thickness";
      parameter Types.Area A_s = (2 * w + 2 * l) * h "elevator case surface area";
      parameter Types.Area A_lift = 0.05 * 0.05 "bucket surface area";
      parameter Types.Area P_lift = 2 * 0.05 + 2 * 0.05 "bucket perimeter";
      parameter Types.Convection h_conv_sw = 180 "solid-wall convection coefficient";
      parameter Types.Convection h_loss = 10 "elevator heat loss coefficient";
      parameter Real dx = h / (n - 2) "discretization size";
      // Particle Variables
      Types.MassFlowRate m_dot_s_in "inlet mass flow rate";
      Types.MassFlowRate m_dot_s_out "outlet mass flow rate";
      Types.Temperature T_s_in "inlet temperature";
      Types.Temperature T_s_out "outlet temperature";
      Types.Temperature T_s[n](each start = 273.15 + 25) "solid temperature distribution";
      Types.Enthalpy h_s[n] "particle enthalpy";
      Types.SpecificHeatCapacity cp_s "particle specific heat";
      Types.Density rho_s "particle density";
      parameter Real phi_s = 0.6 "solid volume fraction";
      Types.Velocity v_s "solid velocity";
      Types.Heat Q_loss "heat loss";
      Types.Temperature T_m[n](each start = 273.15 + 25) "metal temperature distribution";
      Types.SpecificHeatCapacity cp_m "metal specific heat";
      Types.Density rho_m "metal density";
      FallingParticleReceiverSystem.Interfaces.ParticleFlow ParticleOutlet annotation(
        Placement(visible = true, transformation(origin = {0, -100}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-60, 56}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      FallingParticleReceiverSystem.Interfaces.ParticleFlow ParticleInlet annotation(
        Placement(visible = true, transformation(origin = {0, 98}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-60, -58}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
// Connection
      m_dot_s_in = ParticleInlet.m_dot;
      m_dot_s_out = -ParticleOutlet.m_dot;
      T_s_out = ParticleOutlet.T;
      T_s_in = ParticleInlet.T;
// Mass Balance
      m_dot_s_in = m_dot_s_out;
      m_dot_s_in = rho_s * phi_s * v_s * A_lift;
// Energy Balance
      T_s[1] = T_s_in;
      for i in 2:n - 1 loop
        cp_s * rho_s * phi_s * der(T_s[i]) = (-cp_s * rho_s * phi_s * v_s * (T_s[i] - T_s[i - 1]) / dx) + h_conv_sw / P_lift * (T_m[i] - T_s[i]);
      end for;
      T_s[n] = T_s[n - 1];
      T_s[n] = T_s_out;
      T_m[1] = T_m[2];
      for i in 2:n - 1 loop
        cp_m * rho_m * der(T_m[i]) = h_loss / t_m * (T_o - T_m[i]) + P_lift * h / A_s * h_conv_sw / t_m * (T_s[i] - T_m[i]);
      end for;
      T_m[n] = T_m[n - 1];
      Q_loss = sum(h_loss * (T_m[i] - T_o) * dx for i in 2:n - 1) * (A_s / h);
// Properties
      for i in 1:n loop
        h_s[i] = Media.Particle.Enthalpy(T_s[i]);
      end for;
      cp_s = Media.Particle.SpecificHeat(T_s_out);
      rho_s = Media.Particle.Density();
      cp_m = Media.StainlessSteel.SpecificHeat();
      rho_m = Media.StainlessSteel.Density();
      annotation(
        Icon(coordinateSystem(preserveAspectRatio = false)),
        Diagram(coordinateSystem(preserveAspectRatio = false)),
        __OpenModelica_commandLineOptions = "");
    end BucketElevator;

    model BufferVolume
      //Geometric parameters
      parameter Types.Mass m_s_0 = 1 "initial mass inventory";
      parameter Types.Temperature T_s_0 = 20 + 273.15 "initial temperature";
      parameter Types.Volume V_bin = 1 "storage bin volume";
      //Particle variables
      Types.MassFlowRate m_dot_s_in "inlet mass flow rate";
      Types.MassFlowRate m_dot_s_out "outlet mass flow rate";
      Types.Temperature T_s_in "inlet temperature";
      Types.Temperature T_s_out(start = T_s_0) "outlet temperature";
      Types.Enthalpy h_s_in "inlet enthalpy";
      Types.Enthalpy h_s_out "outlet enthalpy";
      Types.Density rho_s "density";
      Types.SpecificHeatCapacity cp_s "specific heat capacity";
      parameter Real phi_s = 0.6 "solids volume fraction";
      Types.Mass m_s(start = m_s_0) "storage bin mass inventory";
      Types.Fraction X_fill "fill fraction";
      FallingParticleReceiverSystem.Interfaces.ParticleFlow ParticleOutlet annotation(
        Placement(visible = true, transformation(origin = {0, -100}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, -100}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      FallingParticleReceiverSystem.Interfaces.ParticleFlow ParticleInlet annotation(
        Placement(visible = true, transformation(origin = {0, 98}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, 100}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
//Connection
      m_dot_s_in = ParticleInlet.m_dot;
      m_dot_s_out = -ParticleOutlet.m_dot;
      T_s_out = ParticleOutlet.T;
      T_s_in = ParticleInlet.T;
//Mass Balance
      der(m_s) = m_dot_s_in - m_dot_s_out;
      X_fill = m_s / (rho_s * phi_s) / V_bin;
      m_dot_s_out = 10 * X_fill ^ 2;
//Energy Balance
      h_s_out * der(m_s) + cp_s * m_s * der(T_s_out) = m_dot_s_in * h_s_in - m_dot_s_out * h_s_out;
//Properties
      h_s_in = Media.ParticleEnthalpy(T_s_in);
      h_s_out = Media.ParticleEnthalpy(T_s_out);
      cp_s = Media.ParticleSpecificHeat();
      rho_s = Media.ParticleDensity();
      annotation(
        Icon(graphics = {Polygon(fillColor = {200, 198, 200}, fillPattern = FillPattern.Solid, points = {{-100, 100}, {-100, -40}, {0, -100}, {100, -40}, {100, 100}, {100, 100}, {-100, 100}}), Polygon(origin = {0, -7}, fillPattern = FillPattern.Solid, points = {{-92, 67}, {0, 93}, {92, 69}, {92, -29}, {0, -85}, {-92, -29}, {-92, 67}})}, coordinateSystem(initialScale = 0.1)));
    end BufferVolume;

    model FallingParticleReceiver1D
      extends Icons.CavityReceiver;
      //Geometry Variables
      parameter Integer n = 100 "number of discretizations";
      parameter Real dx = A_ap ^ 0.5 / (n - 2) "discretization size";
      parameter Types.Temperature T_o = 20 + 273.15;
      parameter Real g = 9.81;
      Types.Length t_c "curtain width";
      parameter Real th_w = 0.05;
      //Particle Distribution Variables
      parameter Types.Temperature T_s_0 = 20 + 273.15 "initial temperature";
      Types.Temperature T_s[n] "solid temperature distribution";
      Types.Velocity v_s[n] "solid velocity";
      Types.Fraction phi_s[n] "solid volume fraction";
      Types.Enthalpy h_s[n] "solid enthalpy";
      Types.SpecificHeatCapacity cp_s[n] "solid specific heat capacity";
      //Particle Variables
      Types.MassFlowRate m_dot_s_in "inlet mass flow rate";
      Types.MassFlowRate m_dot_s_out "outlet mass flow rate";
      Types.Temperature T_s_in "inlet particle temperature";
      Types.Temperature T_s_out "inlet particle temperature";
      Types.Enthalpy h_s_in "inlet particle enthalpy";
      Types.Enthalpy h_s_out "outlet particle enthalpy";
      Types.Density rho_s "particle density";
      //Receiver Wall
      Types.Temperature T_w[n] "backwall temperature";
      Types.Density rho_w "wall density";
      Types.SpecificHeatCapacity cp_w "wall specific heat capacity";
      Types.Conductivity k_w "wall conductivity";
      //Solar Variables
      Types.HeatFlux q_f_solar "solar flux";
      //Types.Heat Q_rerad "reradiation heat loss";
      //Types.Heat Q_in "receiver net input";
      parameter Types.Area A_ap = 1 "aperture area";
      Types.Fraction eta_rec "receiver efficiency";
      parameter Types.Volume V_rec = 1 "receiver volume";
      parameter Types.Fraction sigma_sb = 5.67 * 10 ^ (-8) "stephan boltzman constant";
      parameter Types.Fraction eps_s = 0.88 "emissivity";
      parameter Types.Fraction abs_s = 0.92 "emissivity";
      FallingParticleReceiverSystem.Interfaces.ParticleFlow ParticleInlet annotation(
        Placement(visible = true, transformation(origin = {0, 80}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, 80}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Interfaces.ParticleFlow ParticleOutlet annotation(
        Placement(visible = true, transformation(origin = {0, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealOutput OutletTemperature annotation(
        Placement(visible = true, transformation(origin = {-80, -78}, extent = {{-10, -10}, {10, 10}}, rotation = 180), iconTransformation(origin = {-80, -60}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealInput MassFlow annotation(
        Placement(visible = true, transformation(origin = {-80, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-80, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    initial equation
      for i in 2:n - 1 loop
        T_s[i] = T_s_0;
        v_s[i] = 0.0;
        phi_s[i] = 0.001;
      end for;
    equation
//Connections
      m_dot_s_in = ParticleInlet.m_dot;
      m_dot_s_out = -ParticleOutlet.m_dot;
      T_s_in = ParticleInlet.T;
      T_s_out = ParticleOutlet.T;
      T_s_out = OutletTemperature;
      m_dot_s_in = MassFlow;
//Mass balance
      phi_s[1] = 0.6;
      for i in 2:n - 1 loop
        rho_s * der(phi_s[i]) = -rho_s * (phi_s[i] * v_s[i] - phi_s[i - 1] * v_s[i - 1]) / dx;
      end for;
      phi_s[n] = phi_s[n - 1];
//Momentum balance
      m_dot_s_in = A_ap ^ 0.5 * t_c * rho_s * v_s[1] * phi_s[1];
      v_s[1] = 0.01;
      for i in 2:n - 1 loop
        rho_s * der(phi_s[i] * v_s[i]) = (-rho_s * (phi_s[i] * v_s[i] ^ 2 - phi_s[i - 1] * v_s[i - 1] ^ 2) / dx) + g * rho_s * phi_s[i];
      end for;
      v_s[n] = v_s[n - 1];
      m_dot_s_out = A_ap ^ 0.5 * t_c * rho_s * v_s[n] * phi_s[n];
//Energy balance
      T_s[1] = T_s_in;
      T_w[1] = T_w[2];
      for i in 2:n - 1 loop
        rho_s * cp_s[i] * phi_s[i] * der(T_s[i]) + rho_s * h_s[i] * der(phi_s[i]) = (-rho_s * (phi_s[i] * v_s[i] * h_s[i] - phi_s[i - 1] * v_s[i - 1] * h_s[i - 1]) / dx) + abs_s * q_f_solar / t_c - eps_s * sigma_sb * T_s[i] ^ 4 / t_c - eps_s * sigma_sb * (T_s[i] ^ 4 - T_w[i] ^ 4) / t_c;
        rho_w * cp_w * der(T_w[i]) = k_w * ((T_w[i + 1] - T_w[i]) / dx - (T_w[i] - T_w[i - 1]) / dx) / dx - eps_s * sigma_sb * (T_w[i] ^ 4 - T_s[i] ^ 4) / th_w;
      end for;
      T_s[n] = T_s[n - 1];
      T_s[n] = T_s_out;
      T_w[n] = T_w[n - 1];
//Performance Metrics
      eta_rec * (q_f_solar * A_ap) = m_dot_s_out * h_s_out - m_dot_s_in * h_s_in;
//Properties
      h_s_in = Media.ParticleEnthalpy(T_s_in);
      h_s_out = Media.ParticleEnthalpy(T_s_out);
      for i in 1:n loop
        h_s[i] = Media.ParticleEnthalpy(T_s[i]);
        cp_s[i] = Media.ParticleSpecificHeat();
      end for;
      rho_s = Media.ParticleDensity();
      rho_w = Media.RSLEDensity();
      k_w = Media.RSLEConductivity();
      cp_w = Media.RSLESpecificHeat();
    end FallingParticleReceiver1D;

    model FunnelFlowStorageBin
      extends Icons.StorageBin;
      //Geometric parameters
      parameter Integer n_r = 100;
      parameter Integer n_x = 100;
      parameter Types.Temperature T_o = 20 + 273.15;
      parameter Types.Mass m_s_0 = 1 "initial mass inventory";
      parameter Types.Temperature T_s_0 = 1 "initial temperature";
      parameter Types.Volume V_bin = 1 "storage bin volume";
      parameter Types.Area A_bin = 1 "storage bin cross section";
      parameter Types.Length H_bin = V_bin / A_bin "storage bin height";
      parameter Types.Length R_bin = sqrt(A_bin / Modelica.Constants.pi) "storage bin height";
      parameter Real dx = H_bin / (n_x - 2);
      parameter Real dr = R_bin / (n_r - 2);
      parameter Types.Convection h_loss = 10;
      parameter Types.Area A_s = 2 * A_bin + H_bin * (2 * sqrt(A_bin / 3.14) * 3.14) "storage bin surface area";
      //Particle variables
      Types.MassFlowRate m_dot_s_in "inlet mass flow rate";
      Types.MassFlowRate m_dot_s_out "outlet mass flow rate";
      Types.Temperature T_s_in "inlet temperature";
      Types.Temperature T_s_out(start = T_s_0) "outlet temperature";
      Types.Enthalpy h_s_in "inlet enthalpy";
      Types.Enthalpy h_s_out "outlet enthalpy";
      Types.Density rho_s "density";
      Types.SpecificHeatCapacity cp_s "specific heat capacity";
      parameter Real phi_s = 0.6 "solids volume fraction";
      Types.Mass m_s(start = m_s_0) "storage bin mass inventory";
      Types.Fraction X_fill "fill fraction";
      Types.Heat Q_loss "heat loss";
      FallingParticleReceiverSystem.Interfaces.ParticleFlow ParticleOutlet annotation(
        Placement(visible = true, transformation(origin = {0, -100}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, -100}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      FallingParticleReceiverSystem.Interfaces.ParticleFlow ParticleInlet annotation(
        Placement(visible = true, transformation(origin = {0, 98}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, 100}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
//Connection
      m_dot_s_in = ParticleInlet.m_dot;
      m_dot_s_out = -ParticleOutlet.m_dot;
      T_s_out = ParticleOutlet.T;
      T_s_in = ParticleInlet.T;
//Mass Balance
      der(m_s) = m_dot_s_in - m_dot_s_out;
      X_fill = m_s / (rho_s * phi_s) / V_bin;
//Energy Balance
      h_s_out * der(m_s) + cp_s * m_s * der(T_s_out) = m_dot_s_in * h_s_in - m_dot_s_out * h_s_out - h_loss * A_s * (T_s_out - T_o);
      Q_loss = h_loss * A_s * (T_s_out - T_o);
//Properties
      h_s_in = Media.ParticleEnthalpy(T_s_in);
      h_s_out = Media.ParticleEnthalpy(T_s_out);
      cp_s = Media.ParticleSpecificHeat();
      rho_s = Media.ParticleDensity();
    end FunnelFlowStorageBin;

    model HeliostatField
      extends Icons.Heliostat;
      parameter Integer n = 150 "number of heliostats";
      parameter Types.Length L = 5;
      parameter Types.Length H = 5;
      parameter Types.Area A_ref = L * H;
      parameter Types.Fraction eta_op = 0.6;
      Types.HeatFlux DNI;
      Types.Heat Q_rec;
      Interfaces.Heat Irradiance annotation(
        Placement(visible = true, transformation(origin = {-96, 14}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Interfaces.Heat Insolation annotation(
        Placement(visible = true, transformation(origin = {0, 98}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, 100}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      Insolation.Q = DNI;
      Q_rec = n * A_ref * DNI * eta_op;
      Q_rec = Irradiance.Q;
    end HeliostatField;


    model WeatherData
      extends Icons.Weather;
      Types.HeatFlux DNI;
      Types.Temperature T_amb;
      Interfaces.Heat Insolation annotation(
        Placement(visible = true, transformation(origin = {-48, -96}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-50, -100}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      Insolation.Q = DNI;
//DNI = 950;
      T_amb = 20 + 273.15;
    end WeatherData;

    model DecomposeTime
      extends Icons.Clock;
      Types.Day day "time in days";
      Types.Hour hour "time in hours";
      Types.Day day_number "number of days since start of simulation";
      Types.Hour hour_day "hour since start of new day";
      Types.Second second_hour "seconds since start of new hour";
    algorithm
// Days and Hours
      day := time / (24 * 3600);
      hour := time / 3600;
      day_number := integer(time / (24 * 3600));
      hour_day := integer((time - day_number * 24 * 3600) / 3600);
      second_hour := integer(time - hour_day * 3600 - day_number * 24 * 3600);
    end DecomposeTime;
  end Components;

  package Icons
    extends Library;

    model HeatExchanger
      annotation(
        Icon(graphics = {Ellipse(fillColor = {255, 85, 0}, fillPattern = FillPattern.Sphere, lineThickness = 0.5, extent = {{-100, 100}, {100, -100}}, endAngle = 360), Line(points = {{0, 100}, {0, -100}}, thickness = 0.5), Line(origin = {0, -0.69}, points = {{-100, 0.694251}, {-60, 0.694251}, {-40, 40.6943}, {0, -39.3057}, {40, 40.6943}, {60, 0.694251}, {100, 0.694251}}, thickness = 0.5)}, coordinateSystem(initialScale = 0.1)));
    end HeatExchanger;

    model PropertyPackage
      annotation(
        Icon(graphics = {Rectangle(lineColor = {108, 108, 108}, fillColor = {189, 189, 189}, fillPattern = FillPattern.HorizontalCylinder, extent = {{-100, 100}, {100, -100}}, radius = 25), Ellipse(fillColor = {130, 130, 130}, fillPattern = FillPattern.Sphere, extent = {{-20, 20}, {20, -20}}, endAngle = 360)}));
    end PropertyPackage;

    model Simulation
      annotation(
        Icon(graphics = {Ellipse(lineColor = {0, 85, 0}, fillColor = {219, 229, 219}, fillPattern = FillPattern.Solid, lineThickness = 1, extent = {{-100, 100}, {100, -100}}, endAngle = 360), Polygon(origin = {6, 0}, fillColor = {0, 85, 0}, fillPattern = FillPattern.Solid, points = {{-40, 40}, {-40, -40}, {40, 0}, {-40, 40}})}));
    end Simulation;

    model SolidRes
      annotation(
        Icon(graphics = {Polygon(fillColor = {85, 85, 85}, fillPattern = FillPattern.HorizontalCylinder, points = {{-100, 60}, {-100, -60}, {20, -60}, {100, 0}, {20, 60}, {-100, 60}})}));
    end SolidRes;

    model GasRes
      annotation(
        Icon(graphics = {Polygon(lineColor = {0, 52, 0}, fillColor = {0, 85, 0}, fillPattern = FillPattern.HorizontalCylinder, points = {{-100, 60}, {20, 60}, {100, 0}, {20, -60}, {-100, -60}, {-100, 60}})}));
    end GasRes;

    model Library
      annotation(
        Icon(graphics = {Rectangle(fillColor = {175, 175, 175}, fillPattern = FillPattern.Sphere, extent = {{-100, 100}, {100, -100}}, radius = 25)}));
    end Library;

    model Function
      annotation(
        Diagram,
        Icon(graphics = {Ellipse(fillColor = {255, 255, 127}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}}, endAngle = 360)}),
        __OpenModelica_commandLineOptions = "");
      annotation(
        Diagram,
        Icon(graphics = {Ellipse(fillColor = {255, 255, 127}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}}, endAngle = 360)}),
        __OpenModelica_commandLineOptions = "");
    end Function;

    model Units
      annotation(
        Icon(graphics = {Ellipse(origin = {-1, 0}, fillColor = {0, 0, 255}, fillPattern = FillPattern.Sphere, extent = {{-99, 100}, {99, -100}}, endAngle = 360), Text(fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}}, textString = "SI")}));
    end Units;

    model StorageBin
      annotation(
        Icon(graphics = {Polygon(fillColor = {207, 207, 155}, fillPattern = FillPattern.Solid, points = {{-100, -100}, {-100, 60}, {-40, 100}, {40, 100}, {100, 60}, {100, -100}, {100, -100}, {-100, -100}}), Rectangle(fillPattern = FillPattern.Sphere, extent = {{-80, -20}, {80, -80}}, radius = 2), Rectangle(origin = {0, 90}, fillColor = {89, 89, 89}, fillPattern = FillPattern.Solid, extent = {{-6, 10}, {6, -20}}), Polygon(fillColor = {89, 89, 89}, fillPattern = FillPattern.Solid, points = {{-4, -80}, {-80, -40}, {-80, 40}, {0, 80}, {80, 40}, {80, -38}, {4, -80}, {-4, -80}}), Rectangle(origin = {0, -90}, fillColor = {89, 89, 89}, fillPattern = FillPattern.Solid, extent = {{-4, 10}, {4, -10}})}, coordinateSystem(initialScale = 0.1)));
    end StorageBin;

    model CavityReceiver
      annotation(
        Icon(graphics = {Rectangle(fillColor = {206, 206, 206}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}}, radius = 25), Rectangle(fillColor = {103, 103, 103}, fillPattern = FillPattern.Sphere, extent = {{-80, 80}, {80, -80}}, radius = 10)}, coordinateSystem(initialScale = 0.1)));
    end CavityReceiver;

    model MassFlowHopper
      annotation(
        Icon(graphics = {Polygon(fillColor = {207, 207, 155}, fillPattern = FillPattern.Solid, points = {{-20, -100}, {-100, 80}, {-40, 100}, {40, 100}, {100, 80}, {20, -100}, {0, -100}, {-20, -100}}), Rectangle(origin = {0, 90}, fillColor = {89, 89, 89}, fillPattern = FillPattern.Solid, extent = {{-6, 10}, {6, -20}}), Polygon(fillColor = {89, 89, 89}, fillPattern = FillPattern.Solid, points = {{-4, -80}, {-4, -80}, {-78, 60}, {0, 88}, {78, 60}, {4, -80}, {4, -80}, {-4, -80}}), Rectangle(origin = {0, -90}, fillColor = {89, 89, 89}, fillPattern = FillPattern.Solid, extent = {{-4, 10}, {4, -10}})}, coordinateSystem(initialScale = 0.1)));
    end MassFlowHopper;

    model FunnelFlowHopper
    end FunnelFlowHopper;

    model InsulatedDuct
      annotation(
        Icon(graphics = {Rectangle(fillColor = {207, 207, 155}, fillPattern = FillPattern.Solid, extent = {{-80, 40}, {80, -40}}, radius = 10), Rectangle(fillColor = {122, 122, 122}, fillPattern = FillPattern.HorizontalCylinder, extent = {{100, 25}, {-100, -25}})}, coordinateSystem(initialScale = 0.1)));
    end InsulatedDuct;

    model Heater
      annotation(
        Icon(graphics = {Ellipse(origin = {0, 60}, fillColor = {255, 0, 0}, fillPattern = FillPattern.Sphere, lineThickness = 0.5, extent = {{-60, 40}, {60, -40}}, endAngle = 360), Ellipse(origin = {0, -60}, fillColor = {255, 0, 0}, fillPattern = FillPattern.Sphere, lineThickness = 0.5, extent = {{-60, 40}, {60, -40}}, endAngle = 360), Rectangle(lineColor = {130, 0, 0}, fillColor = {255, 0, 0}, pattern = LinePattern.None, fillPattern = FillPattern.VerticalCylinder, borderPattern = BorderPattern.Raised, extent = {{-60, 60}, {60, -60}}), Rectangle(fillColor = {255, 0, 0}, lineThickness = 0.5, borderPattern = BorderPattern.Raised, extent = {{-60, 60}, {60, -60}}), Line(origin = {-60, 0}, points = {{0, 60}, {0, -60}}, thickness = 0.5), Line(origin = {60, 0.843373}, points = {{0, 60}, {0, -60}}, thickness = 0.5)}, coordinateSystem(initialScale = 0.1)));
    end Heater;

    model BucketElevator
      annotation(
        Icon(graphics = {Rectangle(rotation = -90, fillColor = {207, 207, 155}, fillPattern = FillPattern.Solid, extent = {{-100, 38}, {100, -40}}, radius = 10), Rectangle(rotation = -90, fillColor = {122, 122, 122}, fillPattern = FillPattern.HorizontalCylinder, extent = {{100, 23}, {-94, -25}}), Polygon(origin = {-44, 74}, fillColor = {122, 122, 122}, fillPattern = FillPattern.Sphere, points = {{20, 20}, {-20, -20}, {20, -20}, {20, 20}}), Polygon(origin = {-44, -78}, rotation = 90, fillColor = {122, 122, 122}, fillPattern = FillPattern.Sphere, points = {{20, 20}, {-20, -20}, {20, -20}, {20, 20}})}, coordinateSystem(preserveAspectRatio = false)),
        Diagram(coordinateSystem(preserveAspectRatio = false)),
        __OpenModelica_commandLineOptions = "");
    end BucketElevator;

    model Tee
      annotation(
        Icon(graphics = {Rectangle(origin = {0, 0}, fillColor = {122, 122, 122}, fillPattern = FillPattern.VerticalCylinder, extent = {{-25, 0}, {25, -100}}), Rectangle(fillColor = {122, 122, 122}, fillPattern = FillPattern.HorizontalCylinder, extent = {{100, 25}, {-100, -25}})}, coordinateSystem(initialScale = 0.1)));
    end Tee;

    model Weather
      annotation(
        Icon(graphics = {Ellipse(origin = {-42, 65}, fillColor = {255, 255, 127}, fillPattern = FillPattern.Solid, extent = {{-58, 35}, {42, -65}}, endAngle = 360), Rectangle(origin = {48, 16}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-8, 64}, {12, -56}}), Rectangle(origin = {48, -44}, fillColor = {255, 0, 0}, fillPattern = FillPattern.Solid, extent = {{-8, 64}, {12, -16}}), Ellipse(origin = {54, -66}, lineColor = {255, 0, 0}, fillColor = {255, 116, 116}, fillPattern = FillPattern.Sphere, extent = {{-24, 16}, {16, -24}}, endAngle = 360), Ellipse(origin = {54, -66}, fillColor = {255, 255, 255}, extent = {{-24, 16}, {16, -24}}, endAngle = 360), Ellipse(origin = {50, 80}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-10, 10}, {10, -10}}, endAngle = 180)}));
    end Weather;

    model Heliostat
      annotation(
        Icon(graphics = {Ellipse(origin = {-39, -32}, fillPattern = FillPattern.Solid, extent = {{-9, 8}, {9, -8}}, endAngle = 360), Rectangle(origin = {-15, -34}, fillPattern = FillPattern.Solid, extent = {{-5, 20}, {11, -56}}), Polygon(origin = {0, 11}, fillColor = {207, 207, 207}, fillPattern = FillPattern.Solid, points = {{-80, 9}, {0, -89}, {80, -11}, {0, 69}, {0, 69}, {-80, 9}})}, coordinateSystem(initialScale = 0.1)));
    end Heliostat;

    model Clock
      annotation(
        Icon(graphics = {Ellipse(fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, lineThickness = 1, extent = {{-100, 100}, {100, -100}}, endAngle = 360), Line(origin = {0, 40}, points = {{0, -40}, {0, 32}}, thickness = 1, arrow = {Arrow.None, Arrow.Filled}, arrowSize = 15), Line(origin = {20, 20}, points = {{-20, -20}, {20, 20}}, thickness = 1, arrow = {Arrow.None, Arrow.Filled}, arrowSize = 15)}),
        Diagram,
        __OpenModelica_commandLineOptions = "");
    end Clock;
  end Icons;

  package Models
  extends Icons.Library;

    model G3P3_System
      extends Icons.Simulation;
      //Types.Fraction DNI;
      Types.Fraction HX_op;
      Types.Fraction HX_m_dot;
      Types.Fraction Heater_m_dot;
      Types.Heat Q_loss_total;
      Types.Heat Q_loss_components;
      Types.Heat Q_loss_ducts;
      Types.Heat E_loss(start = 0.0);
      FallingParticleReceiverSystem.Components.FallingParticleReceiver Receiver(A_ap = 1.2 * 1.2, T_s_out(fixed = true), V_rec = 0.01) annotation(
        Placement(visible = true, transformation(origin = {10, 190}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      FallingParticleReceiverSystem.Components.HeliostatField Field annotation(
        Placement(visible = true, transformation(origin = {-42, 152}, extent = {{-10, 10}, {10, -10}}, rotation = 180)));
      FallingParticleReceiverSystem.Components.FreeFallDuct ReceiverDownComer(T_m(fixed = true), h_loss = 0.6, length = 10, perimeter = 0.785, thickness = 0.003) annotation(
        Placement(visible = true, transformation(origin = {10, 164}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      FallingParticleReceiverSystem.Components.LumpedStorageBin HotBin(A_bin = 3.1415 * 2.5 * 2.5, T_s_0 = 1048.15, T_s_out(fixed = true), V_bin = 60, h_loss = 0.25, m_s(fixed = true), m_s_0 = 1000) annotation(
        Placement(visible = true, transformation(origin = {10, 48}, extent = {{-30, -30}, {30, 30}}, rotation = 0)));
      FallingParticleReceiverSystem.Components.WedgeMassFlowHopper ReceiverInventoryHopper(n = 10) annotation(
        Placement(visible = true, transformation(origin = {10, 220}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      FallingParticleReceiverSystem.Components.LumpedStorageBin ColdBin(A_bin = 3.1415 * 2.5 * 2.5, T_s_0 = 888.15, T_s_out(fixed = true), V_bin = 60, h_loss = 0.25, m_s(fixed = true), m_s_0 = 119000) annotation(
        Placement(visible = true, transformation(origin = {10, -130}, extent = {{-30, -30}, {30, 30}}, rotation = 0)));
      FallingParticleReceiverSystem.Components.ParticleHeatExchanger ParticleHeatExchanger(H = 1.5, N_plate = 33, T_CO2(each fixed = false), T_m(each fixed = false), T_s(each fixed = false), W = 0.6, h_conv_CO2 = 3000, h_conv_sw = 450, hc_s = 0.003, n = 20) annotation(
        Placement(visible = true, transformation(origin = {-80, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      FallingParticleReceiverSystem.Components.BucketElevator Elevator(n = 10) annotation(
        Placement(visible = true, transformation(origin = {130, 123}, extent = {{-18, -23}, {18, 23}}, rotation = 0)));
      FallingParticleReceiverSystem.Components.ElectricalHeater ParticleHeater annotation(
        Placement(visible = true, transformation(origin = {-80, 92}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      FallingParticleReceiverSystem.Components.LumpedStorageBin IntermediateStorage(A_bin = 1 * 1, T_s_0 = 888.15, T_s_out(fixed = true), V_bin = 6, h_loss = 0.1, m_s(fixed = true), m_s_0 = 1) annotation(
        Placement(visible = true, transformation(origin = {60, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      FallingParticleReceiverSystem.Components.SolidSplitter HotBinDiverter annotation(
        Placement(visible = true, transformation(origin = {10, 130}, extent = {{-10, 10}, {10, -10}}, rotation = -90)));
      FallingParticleReceiverSystem.Components.SolidJunction HeatExchangerJunction annotation(
        Placement(visible = true, transformation(origin = {-80, 10}, extent = {{-10, 10}, {10, -10}}, rotation = -90)));
      FallingParticleReceiverSystem.Components.FreeFallDuct HotBinBypass(T_m(fixed = true), h_loss = 0.6, length = 10, perimeter = 0.785, thickness = 0.003) annotation(
        Placement(visible = true, transformation(origin = {60, 90}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      FallingParticleReceiverSystem.Components.SolidSplitter ReceiverDiverter annotation(
        Placement(visible = true, transformation(origin = {10, 250}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      FallingParticleReceiverSystem.Components.FreeFallDuct ReceiverBypass(T_m(fixed = true), h_loss = 0.6, length = 10, perimeter = 0.785, thickness = 0.003) annotation(
        Placement(visible = true, transformation(origin = {-80, 130}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      FallingParticleReceiverSystem.Components.FreeFallDuct HeaterDischarge(T_m(fixed = true), h_loss = 0.6, length = 10, perimeter = 0.785, thickness = 0.003) annotation(
        Placement(visible = true, transformation(origin = {-80, 50}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      FallingParticleReceiverSystem.Components.WedgeMassFlowHopper HeatExchangerHopper(n = 10) annotation(
        Placement(visible = true, transformation(origin = {-80, -58}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      FallingParticleReceiverSystem.Components.FreeFallDuct HeatExchangerDownComer(T_m(fixed = true), h_loss = 0.6, length = 10, perimeter = 0.785, thickness = 0.003) annotation(
        Placement(visible = true, transformation(origin = {-40, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      FallingParticleReceiverSystem.Components.SolidJunction ColdBinJunction annotation(
        Placement(visible = true, transformation(origin = {10, -80}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      FallingParticleReceiverSystem.Components.SolidSplitter IntermediateStorageDiverter annotation(
        Placement(visible = true, transformation(origin = {60, -60}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      FallingParticleReceiverSystem.Components.SolidJunction BucketElevatorJunction annotation(
        Placement(visible = true, transformation(origin = {60, -180}, extent = {{-10, 10}, {10, -10}}, rotation = 0)));
      FallingParticleReceiverSystem.Components.FreeFallDuct ColdBinBypass(T_m(fixed = true), h_loss = 0.6, length = 10, perimeter = 0.785, thickness = 0.003) annotation(
        Placement(visible = true, transformation(origin = {60, -130}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      FallingParticleReceiverSystem.Components.FreeFallDuct ColdBinDischarge(T_m(fixed = true), h_loss = 0.6, length = 10, perimeter = 0.785, thickness = 0.003) annotation(
        Placement(visible = true, transformation(origin = {30, -180}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      FallingParticleReceiverSystem.Components.FreeFallDuct HotBinInlet(T_m(fixed = true), h_loss = 0.6, length = 10, perimeter = 0.785, thickness = 0.003) annotation(
        Placement(visible = true, transformation(origin = {10, 100}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      FallingParticleReceiverSystem.Components.FreeFallDuct BucketElevatorDownComer(T_m(fixed = true), h_loss = 0.6, length = 10, perimeter = 0.785, thickness = 0.003) annotation(
        Placement(visible = true, transformation(origin = {50, 280}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      FallingParticleReceiverSystem.SourceSink.CO2Source CO2Inlet annotation(
        Placement(visible = true, transformation(origin = {-110, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      FallingParticleReceiverSystem.SourceSink.CO2Sink CO2Outlet annotation(
        Placement(visible = true, transformation(origin = {-30, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      FallingParticleReceiverSystem.Components.DecomposeTime DecompTime annotation(
        Placement(visible = true, transformation(origin = {-144, 274}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      FallingParticleReceiverSystem.Components.FreeFallDuct HotBinDischarge(T_m(fixed = true), h_loss = 0.6, length = 10, perimeter = 0.785, thickness = 0.003) annotation(
        Placement(visible = true, transformation(origin = {-30, 10}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      FallingParticleReceiverSystem.Components.FreeFallDuct IntermediateStorageDownComer(T_m(fixed = true), h_loss = 0.6, length = 10, perimeter = 0.785, thickness = 0.003) annotation(
        Placement(visible = true, transformation(origin = {60, -30}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      Modelica.Blocks.Sources.Constant Zero(k = 0) annotation(
        Placement(visible = true, transformation(origin = {-170, 210}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Continuous.LimPID ReceiverPID(Ti = 600, controllerType = Modelica.Blocks.Types.SimpleController.PI, initType = Modelica.Blocks.Types.Init.NoInit, k = 2, yMax = 12, yMin = 0.00001, y_start = 10) annotation(
        Placement(visible = true, transformation(origin = {-128, 212}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant ReceiverSetPoint(k = 775 + 273.15) annotation(
        Placement(visible = true, transformation(origin = {-170, 170}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Add add1(k2 = -1) annotation(
        Placement(visible = true, transformation(origin = {-122, 164}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Tables.CombiTable1Ds weatherCombiTable1Ds(tableOnFile=true, tableName="dniData", fileName=Modelica.Utilities.Files.loadResource("C:\Users\Owner\OneDrive\PhD Research\G3P3SystemModeling\GitFiles\Co-Simulation\ABQ_DNI_Lookup.txt"), columns={2}) annotation(
        Placement(visible = true, transformation(origin = {-158, -38}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  Modelica.Blocks.Sources.ContinuousClock continuousClock annotation(
        Placement(visible = true, transformation(origin = {-164, -196}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  Modelica.Blocks.Math.Add ICsubtraction(k2 = -1)  annotation(
        Placement(visible = true, transformation(origin = {-158, -82}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  Modelica.Blocks.Sources.Constant StartTime(k = 0)  annotation(
        Placement(visible = true, transformation(origin = {-102, -110}, extent = {{-10, -10}, {10, 10}}, rotation = 180)));
  Modelica.Blocks.Math.Gain seconds2hours(k = 1 / 3600)  annotation(
        Placement(visible = true, transformation(origin = {-164, -144}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  Components.WeatherData weatherData annotation(
        Placement(visible = true, transformation(origin = {-38, 230}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(add1.y, ReceiverPID.u_m) annotation(
        Line(points = {{-110, 164}, {-128, 164}, {-128, 200}}, color = {0, 0, 127}));
      connect(Receiver.OutletTemperature, add1.u2) annotation(
        Line(points = {{2, 184}, {-98, 184}, {-98, 144}, {-148, 144}, {-148, 158}, {-134, 158}}, color = {0, 0, 127}));
      connect(ReceiverInventoryHopper.ParticleOutlet, Receiver.ParticleInlet) annotation(
        Line(points = {{10, 210}, {10, 198}}));
      connect(Receiver.ParticleOutlet, ReceiverDownComer.ParticleInlet) annotation(
        Line(points = {{10, 182}, {10, 182}, {10, 174}, {10, 174}}));
      connect(ReceiverPID.y, Receiver.MassFlow) annotation(
        Line(points = {{-116, 212}, {2, 212}, {2, 196}, {2, 196}}, color = {0, 0, 127}));
      connect(Zero.y, ReceiverPID.u_s) annotation(
        Line(points = {{-158, 210}, {-142, 210}, {-142, 212}, {-140, 212}}, color = {0, 0, 127}));
      connect(ReceiverSetPoint.y, add1.u1) annotation(
        Line(points = {{-158, 170}, {-134, 170}, {-134, 170}, {-134, 170}}, color = {0, 0, 127}));
      connect(HotBinInlet.ParticleOutlet, HotBin.ParticleInlet) annotation(
        Line(points = {{10, 90}, {10, 78}}));
      connect(HotBin.ParticleOutlet, HotBinDischarge.ParticleInlet) annotation(
        Line(points = {{10, 18}, {10, 10}, {-20, 10}}));
      connect(CO2Inlet.CO2Outlet, ParticleHeatExchanger.CO2Inlet) annotation(
        Line(points = {{-100, -30}, {-90, -30}}));
      connect(ParticleHeatExchanger.CO2Outlet, CO2Outlet.CO2Inlet) annotation(
        Line(points = {{-70, -30}, {-40, -30}, {-40, -30}, {-40, -30}}));
      connect(ParticleHeatExchanger.ParticleOutlet, HeatExchangerHopper.ParticleInlet) annotation(
        Line(points = {{-80, -40}, {-80, -48}}));
      connect(HeatExchangerJunction.ParticleOutlet, ParticleHeatExchanger.ParticleInlet) annotation(
        Line(points = {{-80, -1.49012e-07}, {-80, -20}}));
      connect(IntermediateStorageDownComer.ParticleOutlet, IntermediateStorageDiverter.ParticleInlet) annotation(
        Line(points = {{60, -40}, {60, -40}, {60, -50}, {60, -50}}));
      connect(IntermediateStorage.ParticleOutlet, IntermediateStorageDownComer.ParticleInlet) annotation(
        Line(points = {{60, 0}, {60, 0}, {60, -20}, {60, -20}}));
      connect(HotBinDischarge.ParticleOutlet, HeatExchangerJunction.ParticleInlet2) annotation(
        Line(points = {{-40, 10}, {-70, 10}, {-70, 10}, {-70, 10}}));
      connect(BucketElevatorDownComer.ParticleInlet, Elevator.ParticleOutlet) annotation(
        Line(points = {{60, 280}, {100, 280}, {100, 136}, {120, 136}, {120, 136}}));
      connect(BucketElevatorDownComer.ParticleOutlet, ReceiverDiverter.ParticleInlet) annotation(
        Line(points = {{40, 280}, {10, 280}, {10, 260}, {10, 260}}));
      connect(HotBinDiverter.ParticleOutlet1, HotBinInlet.ParticleInlet) annotation(
        Line(points = {{10, 120}, {10, 120}, {10, 110}, {10, 110}}));
      connect(ColdBin.ParticleOutlet, ColdBinDischarge.ParticleInlet) annotation(
        Line(points = {{10, -160}, {10, -160}, {10, -180}, {20, -180}, {20, -180}}));
      connect(ColdBinDischarge.ParticleOutlet, BucketElevatorJunction.ParticleInlet1) annotation(
        Line(points = {{40, -180}, {50, -180}, {50, -180}, {50, -180}}));
      connect(ColdBinBypass.ParticleOutlet, BucketElevatorJunction.ParticleInlet2) annotation(
        Line(points = {{60, -140}, {60, -140}, {60, -170}, {60, -170}}));
      connect(IntermediateStorageDiverter.ParticleOutlet1, ColdBinBypass.ParticleInlet) annotation(
        Line(points = {{60, -70}, {60, -70}, {60, -120}, {60, -120}}));
      connect(Elevator.ParticleInlet, BucketElevatorJunction.ParticleOutlet) annotation(
        Line(points = {{119.2, 109.66}, {99.2, 109.66}, {99.2, -180.34}, {69.2, -180.34}, {69.2, -180.34}}));
      connect(IntermediateStorageDiverter.ParticleOutlet2, ColdBinJunction.ParticleInlet1) annotation(
        Line(points = {{50, -60}, {10, -60}, {10, -70}, {10, -70}}));
      connect(ColdBinJunction.ParticleOutlet, ColdBin.ParticleInlet) annotation(
        Line(points = {{10, -90}, {10, -90}, {10, -100}, {10, -100}}));
      connect(HeatExchangerDownComer.ParticleOutlet, ColdBinJunction.ParticleInlet2) annotation(
        Line(points = {{-30, -80}, {0, -80}, {0, -80}, {0, -80}}));
      connect(HeatExchangerHopper.ParticleOutlet, HeatExchangerDownComer.ParticleInlet) annotation(
        Line(points = {{-80, -68}, {-80, -80}, {-50, -80}}));
      connect(HeaterDischarge.ParticleOutlet, HeatExchangerJunction.ParticleInlet1) annotation(
        Line(points = {{-80, 40}, {-80, 20}}));
      connect(ParticleHeater.ParticleOutlet, HeaterDischarge.ParticleInlet) annotation(
        Line(points = {{-80, 82}, {-80, 82}, {-80, 60}, {-80, 60}}));
      connect(ReceiverBypass.ParticleOutlet, ParticleHeater.ParticleInlet) annotation(
        Line(points = {{-80, 120}, {-80, 120}, {-80, 102}, {-80, 102}}));
      connect(ReceiverDiverter.ParticleOutlet2, ReceiverBypass.ParticleInlet) annotation(
        Line(points = {{-1.49012e-07, 250}, {-80, 250}, {-80, 140}, {-80, 140}}));
      connect(ReceiverDiverter.ParticleOutlet1, ReceiverInventoryHopper.ParticleInlet) annotation(
        Line(points = {{10, 240}, {10, 240}, {10, 230}, {10, 230}}));
      connect(HotBinBypass.ParticleOutlet, IntermediateStorage.ParticleInlet) annotation(
        Line(points = {{60, 80}, {60, 20}}));
      connect(HotBinDiverter.ParticleOutlet2, HotBinBypass.ParticleInlet) annotation(
        Line(points = {{20, 130}, {60, 130}, {60, 100}}));
      connect(ReceiverDownComer.ParticleOutlet, HotBinDiverter.ParticleInlet) annotation(
        Line(points = {{10, 154}, {10, 154}, {10, 140}, {10, 140}}));
  connect(StartTime.y, ICsubtraction.u2) annotation(
        Line(points = {{-113, -110}, {-152, -110}, {-152, -94}}, color = {0, 0, 127}));
  connect(ICsubtraction.y, weatherCombiTable1Ds.u) annotation(
        Line(points = {{-158, -71}, {-158, -50}}, color = {0, 0, 127}));
  connect(continuousClock.y, seconds2hours.u) annotation(
        Line(points = {{-164, -184}, {-164, -156}}, color = {0, 0, 127}));
  connect(seconds2hours.y, ICsubtraction.u1) annotation(
        Line(points = {{-164, -132}, {-164, -94}}, color = {0, 0, 127}));
  connect(weatherData.Insolation, Field.Insolation) annotation(
        Line(points = {{-42, 220}, {-42, 162}}));
    equation
      weatherData.DNI = weatherCombiTable1Ds.y[1] + 1E-10;
//  DNI = weatherData.DNI;
// Heat Exchanger Flow Rate
      CO2Inlet.T = 565 + 273.15;
      CO2Inlet.m_dot = 5;
// Receiver Flow Rate
      if ColdBin.m_s > 1000 then
        Receiver.Q_solar = Field.Q_rec;
      else
        Receiver.Q_solar = 0.0001 * Field.Q_rec;
      end if;
// Receiver Diverter Valve Control
      if HotBinDiverter.T_s_in > 765 + 273.15 and HotBinDiverter.T_s_in < 800 + 273.15 then
        HotBinDiverter.m_dot_s_out_1 = HotBinDiverter.m_dot_s_in;
      else
        HotBinDiverter.m_dot_s_out_2 = HotBinDiverter.m_dot_s_in;
      end if;
// Intermediate Storage Diverter Valve Control
      IntermediateStorageDiverter.m_dot_s_out_2 = 0;
// Heat Exchanger Flow Rate
      HotBin.m_dot_s_out = HX_op * HX_m_dot;
      if DecompTime.hour_day > 2.0 and DecompTime.hour_day < 8.0 then
        HX_op = 1.0;
      else
        HX_op = 1E-10;
      end if;
// Electrical Heater Flow Rate
      if HotBin.m_dot_s_out < 0.01 then
        ParticleHeater.m_dot_s_out = Heater_m_dot;
      else
        ParticleHeater.m_dot_s_out = 0.0;
      end if;
      ParticleHeater.T_s_out = 615 + 273.15;
// Total Heat Loss
      Q_loss_components = ColdBin.Q_loss + HotBin.Q_loss + IntermediateStorage.Q_loss + Elevator.Q_loss;
      Q_loss_ducts = BucketElevatorDownComer.Q_loss + ReceiverBypass.Q_loss + ReceiverDownComer.Q_loss + HotBinBypass.Q_loss + HotBinDischarge.Q_loss + HeaterDischarge.Q_loss + HeatExchangerDownComer.Q_loss + IntermediateStorageDownComer.Q_loss + ColdBinDischarge.Q_loss;
      Q_loss_total = Q_loss_components + Q_loss_ducts;
      der(E_loss * 3600) = Q_loss_total / 1E6;
// Control Algorithm
    algorithm
// Heat Exchanger Operation
      when initial() then
        HX_m_dot := 1E-10;
        Heater_m_dot := 0.1;
      elsewhen HotBin.m_s > 10000 then
        HX_m_dot := 5.0;
      elsewhen HotBin.m_s < 1000 then
        HX_m_dot := 1E-10;
      end when;
// Receiver Operation
//  when initial() then
//    DNI := 1E-10;
//  elsewhen DecompTime.hour_day > 10.0 then
//    if DecompTime.day_number == 5.0 then
//      DNI := 1E-10;
//    elseif DecompTime.day_number == 3.0 then
//      DNI := 900;
//    elseif DecompTime.day_number == 7.0 then
//      DNI := 850;
//   else
//      DNI := 1000;
//    end if;
//  elsewhen DecompTime.hour_day > 16.0 then
//    DNI := 1E-10;
//  end when;
// Intermediate Storage Bin Discharge Control
      when initial() then
        IntermediateStorage.m_dot_s_out := 0.0;
      elsewhen IntermediateStorage.m_s > 10 then
        IntermediateStorage.m_dot_s_out := 1.0;
      elsewhen IntermediateStorage.m_s < 1.0 then
        IntermediateStorage.m_dot_s_out := 0.0;
      end when;
// Energy Loss Reset
//  when DecompTime.hour_day > 23.999 then
//    reinit(E_loss, 0.0);
//  end when;
      annotation(
        Diagram(coordinateSystem(extent = {{-200, -300}, {200, 300}})),
        Icon,
        __OpenModelica_commandLineOptions = "");
    end G3P3_System;

    model LumpedStorageBin_Demo
      extends Icons.Simulation;
      FallingParticleReceiverSystem.Components.LumpedStorageBin lumpedStorageBin1(A_bin = 20, T_s_0 = 1048.15, T_s_out(fixed = true), V_bin = 60, m_s(fixed = true), m_s_0 = 1000) annotation(
        Placement(visible = true, transformation(origin = {0, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      SourceSink.ParticleSource particleSource1 annotation(
        Placement(visible = true, transformation(origin = {-30, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      FallingParticleReceiverSystem.SourceSink.ParticleSink particleSink1 annotation(
        Placement(visible = true, transformation(origin = {30, -20}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      FallingParticleReceiverSystem.Components.DecomposeTime Time annotation(
        Placement(visible = true, transformation(origin = {-70, 70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(lumpedStorageBin1.ParticleOutlet, particleSink1.ParticleInlet) annotation(
        Line(points = {{0, 0}, {0, -20}, {20, -20}}));
      connect(particleSource1.ParticleOutlet, lumpedStorageBin1.ParticleInlet) annotation(
        Line(points = {{-20, 40}, {0, 40}, {0, 20}, {0, 20}}));
      particleSource1.T = 775 + 273.15;
      if time < 3600 * 6 then
        lumpedStorageBin1.m_dot_s_out = 0;
        lumpedStorageBin1.m_dot_s_in = 5;
      else
        lumpedStorageBin1.m_dot_s_out = 5;
        lumpedStorageBin1.m_dot_s_in = 0;
      end if;
    end LumpedStorageBin_Demo;

    model LumpedParticleReceiver_Demo
      extends Icons.Simulation;
      FallingParticleReceiverSystem.Components.FallingParticleReceiver Receiver annotation(
        Placement(visible = true, transformation(origin = {0, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      FallingParticleReceiverSystem.SourceSink.ParticleSource ReceiverInlet annotation(
        Placement(visible = true, transformation(origin = {-30, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      FallingParticleReceiverSystem.SourceSink.ParticleSink ReceiverOutlet annotation(
        Placement(visible = true, transformation(origin = {30, -40}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    equation
      connect(Receiver.ParticleOutlet, ReceiverOutlet.ParticleInlet) annotation(
        Line(points = {{0, 2}, {0, -40}, {20, -40}}));
      connect(ReceiverInlet.ParticleOutlet, Receiver.ParticleInlet) annotation(
        Line(points = {{-20, 60}, {0, 60}, {0, 18}}));
    equation
      ReceiverInlet.T = 600 + 273.15;
      ReceiverInlet.m_dot = 10;
      Receiver.Q_solar = 1000000;
    end LumpedParticleReceiver_Demo;

    model HeliostatField_Demo
      extends Icons.Simulation;
      FallingParticleReceiverSystem.Components.WeatherData Weather annotation(
        Placement(visible = true, transformation(origin = {0, 70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      FallingParticleReceiverSystem.Components.HeliostatField Field annotation(
        Placement(visible = true, transformation(origin = {0, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(Weather.Insolation, Field.Insolation) annotation(
        Line(points = {{-4, 60}, {0, 60}, {0, 20}, {0, 20}}));
    end HeliostatField_Demo;

    model ElectricalHeater_Demo
      extends Icons.Simulation;
    end ElectricalHeater_Demo;
  end Models;
  annotation(
    uses(Modelica(version = "3.2.2")));
end FallingParticleReceiverSystem;

model exampleTableLookup
  package Medium=Modelica.Media.IdealGases.SingleGases.CO2;
  inner Modelica.Fluid.System system(energyDynamics=Modelica.Fluid.Types.Dynamics.FixedInitial) annotation(
    Placement(visible = true, transformation(origin = {-66, 76}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Tables.CombiTable1Ds combiTable1Ds(tableOnFile=true, tableName="dniData", fileName=Modelica.Utilities.Files.loadResource("C:\Users\Owner\OneDrive\PhD Research\Modelica Model\ABQ_DNI_Lookup.txt"), columns={2}) annotation(
    Placement(visible = true, transformation(origin = {50, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.ContinuousClock continuousClock annotation(
    Placement(visible = true, transformation(origin = {-62, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Gain second2hour(k = 1/3600)  annotation(
    Placement(visible = true, transformation(origin = {10, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
equation
  connect(continuousClock.y, second2hour.u) annotation(
    Line(points = {{-50, 30}, {-2, 30}}, color = {0, 0, 127}));
  connect(second2hour.y, combiTable1Ds.u) annotation(
    Line(points = {{22, 30}, {38, 30}}, color = {0, 0, 127}));
  annotation(
    uses(Modelica(version = "4.0.0")),
    Diagram);
end exampleTableLookup;

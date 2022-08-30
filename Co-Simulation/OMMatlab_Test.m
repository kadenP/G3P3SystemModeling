% verifies that OMMatlab and Modelica are working properly

import OMMatlab.*
omc = OMMatlab();

% identify version of Modelica being used
ModelicaVer = omc.sendExpression("getVersion()")

% simulate simple sine function using the Modelica engine
omc.sendExpression("loadModel(Modelica)");
omc.sendExpression("model a Real s; equation s=sin(10*time); end a;");
omc.sendExpression("simulate(a)");
omc.sendExpression("plot(s)");


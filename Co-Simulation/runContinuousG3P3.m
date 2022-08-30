% runs the G3P3 model as is (without heat kernel model)

import OMMatlab.*
omc = OMMatlab();

% load G3P3 Modelica Model
omc.ModelicaSystem("FallingParticleReceiverSystem-REV4.mo", ...
    "FallingParticleReceiverSystem.Models.G3P3_System", ["Modelica", "Complex", ...
    "ModelicaServices", "ModelicaReference"]);

% check avaliable model parameters
% omc.getParameters()

% setup and run simulation
% omc.setSimulationOptions(["stopTime=259200", "stepSize=60", ...
%                           "solver=dassl", "tolerance=1e-08"]);
omc.setSimulationOptions(["stopTime=180", "stepSize=60", ...
                          "solver=dassl", "tolerance=1e-08"]);
omc.simulate("G3P3_results.mat");

% get the simulation results
availableResults = omc.getSolutions();   % list of resulting parameters
results = omc.getSolutions(["time", "DecompTime.hour", "HotBin.T_s_out", "HotBin.T_s_in", ...
                           "HotBin.m_s", "HotBin.m_dot_s_out", "HotBin.m_dot_s_in", "HotBin.X_fill", "IntermediateStorage.m_s", ...
                           "IntermediateStorage.m_dot_s_out", "IntermediateStorage.m_dot_s_in", ...
                           "DNI", "Field.DNI", "Weather.DNI"]);
                       
HX_params = omc.getSolutions(["DecompTime.hour", ...
                                   "DNI", ...
                                   "ParticleHeater.Q_in", ...
                                   "ParticleHeater.m_dot_s_out", ...
                                   "ParticleHeatExchanger.m_dot_CO2_in", ...
                                   "ParticleHeatExchanger.m_dot_CO2_out", ...
                                   "ParticleHeatExchanger.m_dot_s_in", ...
                                   "ParticleHeatExchanger.m_dot_s_out"]);  

% get resulting parameters
endParams = omc.getParameters();
fns = fieldnames(endParams);

% % test methodology for setting new parameters with old (new IC)
% newIC = [];
% for i = 1:length(fns)
%     if ~isempty(endParam.(fns{i}))
%         newIC = [newIC, convertCharsToStrings(strcat(fns{i}, '=', endParam.(fns{i})))];
%     end
% end
% 
% 
% omc.setParameters(newIC);
% omc.getParameters()












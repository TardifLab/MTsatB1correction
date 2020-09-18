function [R1app_vfa, Aapp_vfa] = MAMT_model_simVFA(Params)
%% Recreate the MAMT Model from Portnoy and Stanisz (2007)

% Sim the VFA experiment for the calculation of MTsat.
%% Currently assumes 30ms TR and 5 and 20 degree flip angle. 
Params.b1 = 0; % microTesla
Params.numSatPulse =1;
Params.pulseDur = 1/1000; %duration of 1 MT pulse in seconds
Params.GapDur = 0.1/1000; %ms gap between MT pulses in train
Params.TR = 30/1000; % total repetition time = MT pulse train and readout.
Params.WExcDur = 3/1000; % duration of water pulse
Params.numExcitation = 1; % number of readout lines/TR
Params.flipAngle = 5; % excitation flip angle water.

lfa = MAMT_model_2007_5(Params);

% second flip angle
Params.flipAngle = 20; % excitation flip angle water.

hfa = MAMT_model_2007_5(Params);

%% VFA R1 and Aapp calculation
a1 = 5 * pi/180;
a2 = 20 * pi / 180;
TR = Params.TR;
R1app_vfa = 0.5 .* (hfa.*a2./ TR - lfa.*a1./TR) ./ (lfa./(a1) - hfa./(a2));
Aapp_vfa = lfa .* hfa .* (TR .* a2./a1 - TR.* a1./a2) ./ (hfa.* TR .*a2 - lfa.* TR .*a1);




























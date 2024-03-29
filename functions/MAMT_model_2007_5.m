function f_s = MAMT_model_2007_5(Params)
%% Similar to the MAMT Model from Portnoy and Stanisz (2007)

% %% Required Pool values
% Params.M0a = 1;
% Params.M0b = 0.01;
% Params.R = 0; % exchange rate;
% Params.b1 = 0; % microTesla
% Params.Raobs = 0.431;
% Params.Rb = 1;
% Params.T2a = 45e-3; 
% Params.T2b = 12e-6; 
% Params.delta = 5000;
% Params.T1D = 1e-3;
% 
% % Sequence Parameters
% Params.pulseTrainLength =1;
% Params.pulseDur = 5; %ms duration of 1 MT pulse
% Params.pulseGapDur = 0.25/1000; %ms gap between MT pulses in train
% Params.TR = 8; % total repetition time = MT pulse train and readout.
% Params.WExcDur = 3/1000; % duration of water pulse
% Params.numExcitation = 1; % number of readout lines/TR
% Params.flipAngle = 90; % excitation flip angle water.
% Params.WExcw1 = Params.flipAngle/(360*Params.WExcDur); % not needed anymore.
% Params.lineshape = 'gaussian'; % 'gaussian' or 'superLor';
% Params.SatPulseShape = 'hanning'; % options: 'hanning', 'gaussian', 'square'

stepSize = 50e-6; % 50 microseconds
Params.stepSize = stepSize;

%% Initialize
ref_delta = Params.delta;

if isempty(Params.Ra) % allow you to specify either Ra or Raobs
    Params.Ra = Params.Raobs - ((Params.R * Params.M0b * (Params.Rb - Params.Raobs)) / (Params.Rb - Params.Raobs + Params.R));
    if isnan(Params.Ra)
        Params.Ra = 1;
    end
end

Bpr0 = 0; % define for consistency

% M is a 2D Magnetization Vector containing Mza, Mzb
M0 = [Params.M0a, Params.M0b, Bpr0]';
I = eye(3); % identity matrix
        
B = [Params.Ra*Params.M0a, Params.Rb*Params.M0b, Bpr0]';

%% Calculate timing variables
loops = 600;

PulseDur = ceil(Params.pulseDur/stepSize); % iteration number
%GapDur = Params.pulseGapDur/stepSize; % iteration number
WExcDur = ceil(Params.WExcDur/stepSize); % iteration number
TRDur = round( Params.TR/stepSize); % iteration number

% enable multi-echo readout
if Params.numExcitation > 1
    echoSpacing = Params.echospacing; % time in seconds
else
    echoSpacing = 0;
end

TR_fill = Params.TR - (Params.numSatPulse)*( Params.pulseDur + Params.pulseGapDur) - Params.numExcitation*( Params.WExcDur + echoSpacing); % time in seconds
%TR_fillDur = TR_fill/stepSize;

%% Calculation
M_t = zeros(3,TRDur*loops +1);
M_t(:,1) = M0;
idx = 2;
prev_val = 0;

% calculate time-vary RF pulses...
[satPulse, sinc_ExcPulse] = MAMT_preparePulses(Params);

for i = 1:loops
    %m %report progress...
    for j = 1:Params.numSatPulse % for each MTsat pulse cycle

        %%For the MT pulse

        for k = 1:PulseDur
            
            if strcmp( Params.freqPattern,'dualAlternate')
                if (mod(j,2) == 0) % alternate every other pulse
                    Params.delta = -ref_delta;
                else
                    Params.delta = ref_delta;
                end
                
            else % single is default;   
                Params.delta = ref_delta; % dualContinuous handled in matrix calc
            end
                        
            Params.w1 = 2*pi *satPulse(k) .* 42.57747892;
            t = stepSize;

            A = calc_RF_matrix_wDipolar2(Params); % update the RF values in matrix
            M_t(:,idx) = expm(A*t) * M_t(:,idx-1) + (expm(A*t) - I)* (A\B); % Update Magnetization.                           
            
            idx = idx +1;    
        end    
   
        %%During the Pulse Gap, the pools relax, Rrfb = 0
        for k = 1 %:GapDur
            Params.w1 = 0;
            Params.delta = 0;
            t = Params.pulseGapDur; % calculate the whole relaxation at once for speed

            A = calc_RF_matrix_wDipolar2(Params); % update the RF values in matrix
            M_t(:,idx) = expm(A*t) * M_t(:,idx-1) + (expm(A*t) - I)* (A\B); % Update Magnetization.                           
            
            if (j == Params.numSatPulse) % && (k == GapDur)
                check_val = M_t(1,idx);

                diff_val = abs(check_val - prev_val)*100;

                if  diff_val < 0.05
                    f_s = M_t(1,idx)*sin(Params.flipAngle *pi/180); % if it has hit Steady state, finish. 
                    return;
                else
                    prev_val = check_val;
                end
            end
            
            idx = idx +1;
                
        end                
    end
    
    ex_count = 0;
    for j = 1: Params.numExcitation
        for k = 1:WExcDur % pool relaxation during the gap
            Params.w1 = sinc_ExcPulse(k) *42.57748; %Params.WExcw1;
            Params.delta = 0;
            t = stepSize;
            
            A = calc_RF_matrix_wDipolar2(Params); % update the RF values in matrix
            M_t(:,idx) = expm(A*t) * M_t(:,idx-1) + (expm(A*t) - I)* (A\B); % Update Magnetization.     
            
            idx = idx +1;   
            if k == WExcDur
                ex_count = ex_count +1;
            end
        
        end
        
        if ex_count < Params.numExcitation
            for k = 1 % pool relaxation during the gap
                Params.w1 = 0;
                Params.delta = 0;
                t = echoSpacing; % calculate the whole relaxation at once for speed
                
                A = calc_RF_matrix_wDipolar2(Params); % update the RF values in matrix
                M_t(:,idx) = expm(A*t) * M_t(:,idx-1) + (expm(A*t) - I)* (A\B); % Update Magnetization.
        
                idx = idx +1;  
            end
        end
        

    end
    
    for k = 1 % :TR_fillDur % pool relaxation during the gap
        Params.w1 = 0;
        Params.delta = 0;
        t = TR_fill; % calculate the whole relaxation at once for speed

        A = calc_RF_matrix_wDipolar2(Params); % update the RF values in matrix
        M_t(:,idx) = expm(A*t) * M_t(:,idx-1) + (expm(A*t) - I)* (A\B); % Update Magnetization.

        idx = idx +1;  
    end

end


%% Use this to trouble shoot and visualize the magnetization vectors
% M_t(:,idx:end) = [];
% 
% figure; plot(M_t(1,:))
% figure; plot(M_t(2,:))
% 



















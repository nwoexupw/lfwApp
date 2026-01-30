
% Date: 9/25/24 - 10/2/24
% Author: Vasu Bhardwaj, Gordon Parker , Yagel Blinkoff  
% This code is to generate control desk sdf file to communicate with dspace
% hardware for Low friction testbed (LFT), in LFT test bed we give current/
% force as input in form of feedforward or controller. 

%% Cleanup
clearvars;
close('all');

%% Open Simulink Model
mdl.name = 'lftModel';  % define model name
open_system(mdl.name);       % open Simulink model

%% Laser Sensor Configuration

z_stroke = 89.3207; % Total stroke length (BDC-TDC) [mm]
z_TDC    = 1.9760;  % Sensor reading at TDC [mm]

z_bias      = z_TDC; % subtracted to the laser data, to achive 0mm at TDC
z_midShift  = z_stroke/2; % subtracted to make 0 mm position at the center of the stroke


%% Define Variant for FeedForward and Controller input type
% In this section we are assigning a constant value to every input type
% signal to use in variant.In Feedforward there are 2 type of input current
% and force 

% Feedforward: Current input 
ZERO = 0; % Zero Current input 
CONS = 1; % Constant current input 
SINE = 2; % Sine wave input 
TRAJ = 5; % Trajectory input
WNOIZ = 3; % white Noise input 

% Feedforward: Constants Force
FCONS = 7; % Constant force input 
FTRAJ = 6; % Trajectory force input 
FWNOIZ = 9;  % White noise input 

% Feedbackward: Controller 
DAMP = 8; % Damping controller 
SPED = 4; % 



%% Variant Selection: USER INPUT 
% Select the feedforward variant and signal shape for simulation
FEEDFWD = CONS;   % user need to update it
FEEDBAK = DAMP;   % type of controller using 


% Feedforward (FEEDFWD) : Current input 
% Choose ZERO for Zero Current input 
% Choose CONS for Constant current input 
% Choose SINE for Sine wave current input 
% Choose TRAJ for Trajectory current input
% Choose WNOIZ for white Noise Current input 

% Feedforward (FEEDFWD): Force input
% Choose FCONS for constant force input 
% Chooe FTRAJ for Trajectory force input 
% Choose FWNOIZ for White noise force input 


ff.draft = 0; % N, as per weight of Buoy % Adjust Draft  


%% Controller Settings
% Velocity feedback
fb.k = -0.1; % N/(m/s)

% Damping coefficient 
dampC = 0.1; %

%% Feedforward Setting

if FEEDFWD == SINE
    % Sine
    ff.a = 0.01; % amplitude, 
    ff.f = 1/60; % frequency, Hz


elseif FEEDFWD == TRAJ
    % Trajectory Input
    % User Input 
    freqVector = 0.17 ;%2.07; %Hz
    NumCycle_ramp = 0.1; % Number of single sine periods for ramp up and ramp down
    NumCycle_ss = 20; % Number of single sine periods for steady state (between ramp up and ramp down) 
    TimeStep = 0.01; % Desired time step in seconds
    % CoolDownPeriod = 20; % Cool Down Period in seconds
    CoolDownPeriod = 0.1;
    
    
    
    % Generate ramped single sine signals 
    [signal, time] = genRampSingleSine(freqVector, NumCycle_ramp, NumCycle_ss, TimeStep, CoolDownPeriod);
    
    % Make trajectory
    trajData.t = time;
    trajData.z = signal; % Convert traj to counts
    
    % Plot the generated signals 
    figure 
    plot(time, signal) 
    grid on 
    xlabel('Time [s]') 
    ylabel('Current Signal')



elseif FEEDFWD == WNOIZ
    % White Noise Settings
    wn = makeWhiteNoise(100,30,20); % sample freq (hz), tEnd (s), bandwidth (hz)
    wn.max= 1;
    
    % load('whiteNoiseData_300s.mat')
    % 
    % wn.y = OUT(1).xt; 
    % wn.t = 0:1/100:length(wn.y)*1/100-1/100;



    figure 
    plot(wn.t, wn.y) 
    grid on 
    xlabel('Time [s]') 
    ylabel('Current Signal')


elseif FEEDFWD == FTRAJ
    % The loaded xData and yData can now be used for interpolation or other operations
    % to simulate the system's behavior based on the lookup table.
    
    % %%Trajectory Input
    % 
    % % 6/12/25 - wtf. We just wasted 20 min trying to understand if this is
    % % used anywhere. We think it is not. Someone fix this!!!
    % % Rampsine
    % fIN.freqVector = 0.17; %Hz
    % fIN.NumCycle_ramp = 1; % Number of single sine periods for ramp up and ramp down
    % fIN.NumCycle_ss = 10; % Number of single sine periods for steady state (between ramp up and ramp down) 
    % fIN.TimeStep = 0.01; % Desired time step in seconds
    % % CoolDownPeriod = 20; % Cool Down Period in seconds
    % fIN.CoolDownPeriod = 1;
    % 
    % 
    % % Generate ramped single sine signals 
    % [fIN.signal, fIN.time] = genRampSingleSine(fIN.freqVector, fIN.NumCycle_ramp, fIN.NumCycle_ss, fIN.TimeStep, fIN.CoolDownPeriod);
    % 
    % % Make trajectory
    % fIN.tnew = fIN.time;
    % fIN.Fnew = fIN.signal; % Convert traj to counts

    % % Plot the generated signals 
    % figure 
    % plot(fIN.tnew, fIN.Fnew) 
    % grid on 
    % xlabel('Time [s]') 
    % ylabel('Force Signal')

    ampN = 1; % We always keep this one at 1
    fHz = .085; % Change this frequency (Hz)
    T = 1/fHz; % Ramp time ( 1 cycle)
   

    % Plot the generated signals 
    figure 
    t = 0:0.001:20; 
    plot(t, ampN * sin(2*pi*fHz*t)) 
    grid on 
    xlabel('Time [s]') 
    ylabel('Force Signal')

elseif FEEDFWD == FWNOIZ
     % White Noise Settings
    % fwn = makeWhiteNoise(100,30,20); % sample freq (hz), tEnd (s), bandwidth (hz)
    % fwn.max= 1;
    
    % Load White Noise Data
    %load("whiteNoiseData_300s.mat")
    % load("pinkNoiseData.mat")
    
    % fwn.expNum = 3; % Exp # 
    % fwn.fs = 100; 
    % fwn.y = OUT(fwn.expNum).xt(1:400e2);   
    % fwn.t = (0:1/fwn.fs:length(fwn.y)/fwn.fs-1/fwn.fs)'; 
    % 
    % figure 
    % plot(fwn.t, fwn.y) 
    % grid on 
    % xlabel('Time [s]') 
    % ylabel('Force Signal')

end 



%% Define Variants

% Variant for Feed fordward 
FF_ZERO_VSS = Simulink.Variant('FEEDFWD == ZERO');
FF_CONS_VSS = Simulink.Variant('FEEDFWD == CONS');
FF_SINE_VSS = Simulink.Variant('FEEDFWD == SINE');
FF_NOIZ_VSS = Simulink.Variant('FEEDFWD == WNOIZ');
FF_TRAJ_VSS = Simulink.Variant('FEEDFWD == TRAJ');

FF_FTRAJ_VSS = Simulink.Variant('FEEDFWD == FTRAJ');
FF_FCONS_VSS = Simulink.Variant('FEEDFWD == FCONS');
FF_FWN_VSS = Simulink.Variant('FEEDFWD == FWNOIZ');


% Variant for Controller's 
FB_ZERO_VSS = Simulink.Variant('FEEDBAK == ZERO');
FB_SPED_VSS = Simulink.Variant('FEEDBAK == SPED');
FB_Damp_VSS = Simulink.Variant('FEEDBAK == DAMP');


%% Model Settings and Input 

md.dt = 0.001; % time step, s

% saturation input 
Current_Limit.max = 4; % saturation input, amp
fIN.Forcemax = 30; % saturation force cmd, N

% Note: needs adjustment if using laser frame for position.
fIN.POSmax = z_stroke/2; % saturation position cmd,mm 
fIN.POSmin = -z_stroke/2 ; % saturation position cmd,mm


% LCAM Settings
% amplifier current sense gain (LCAM pin 7), A/V. This is unchangeable. It
% it in the LCAM manual, e.g. if pin 7 shows 1 v, 2 A is flowing.
lcam.v2i = 2.0; 

% amplifier gain set using RV1 (labeled Signal Gain on the LCAM chassis).
% Note: make sure that the Loop Gain (the funny looking pot) is all the
% way CW. This means that the LCAM is ON. Also, make sure that the current
% output is 0 when 0 volts is commanded to Signal+ and Signal-. If it's
% not, then use RV2 (Balance) to zero it out.
%lcam.gain = 0.44; % amplifier gain, A/V 
lcam.gain = 1;  % amplifier gain, A/V 

% Load Cell Settings
% Volts to Newton (Load cell)
lc.k=1.1314;  % gain, g/volt
lc.b=-0.0127; % bias, g
lc.s=9.81;    % scale, GRAVITY

% Laser Displacement Sensor Settings
lz.o = 0;  % distance from laser read head to target when fully UP

RXb = 4;   % number of bytes to recieve
TXb = 2;   % number of dat bytes to send [2 cmd, 1 space, 5 argument]                   
cts = 1;   % clear to send signal

% AR100 (Laser) sensor parameters
ar.range = 100;
ar.domain = 16384;
ar.offset = 93.3; % offset (mm) from AR sensor measurement start to bottom of plate

% Load and Utilize Lookup Table Data (Force as input) 
% Load the pre-generated lookup table data from the MAT file
load("VCC_04_21_2025.mat")

%load('ForceVsPosition_LT_9_9_24.mat');
%load('ForceVsPosition_LT_3_17_25.mat'); % New lookup table
%test_xData = -1* xData; % Switch X axis for better representaion 
% xData = flip(xData); 
% yData = flip(yData); 

% Load the pre-generated lookup table data from the MAT file
%load('LTMegforce9_9_midpointZero.mat');
%load('LookupT_MagF_LT_3_17_25.mat'); % New lookup table
%Mag_xData = -1* Mag_xData; 
% Mag_xData = flip(Mag_xData); 
% Mag_yData = flip(Mag_yData);


%% Helper Function(s)
function wn = makeWhiteNoise(fs,tEnd,bw)
% MAKEWHITENOISE create a a white noise structure. 

  rng default

  % Parameters that determine the sequence length, etc.
  %fs   = 100; % sample rate, e.g. dSPACE update, sample/sec
  %tEnd =  30; % length of data sequence, sec
  %bw   =  10; % data sequence cutoff frequency, hz

  % Create samples
  x  = randn(fs*tEnd,1);   % full spectrum, white noise data
  fx = fft(x);             % frequency domain
  N  = length(fx);         % data sequence length

  % Implement cutoff frequency as percent of frequency bins
  p = bw/fs; 

  % Frequency bin length
  P = round(p*N); 

  % Chop the fft to exclude frequencies above the cutoff
  fx(P+1:end) = 0; 

  % Convert back to time domain
  wn.y = real(ifft(fx));
  wn.t = ((0:length(wn.y)-1)')/fs;
  %wn= timeseries(y, t);
end

%% Sweep Function 


% Create valid Simulink parameters

freq_up = 0.1 : 0.1 : 3;        % Establish increasing frequency range 
freq_down = 3.0 : -0.1 : 0.1;   % Establish decreasing frequency range
freq = [freq_up,freq_down] ;


cyc = zeros(1, length(freq_down));
cyc(:) = 25;                    % Establish vector of cycles for all frequencies

[cycles,freqs] = createSimParams(cyc,freq_down);


% DONE !!

% Helper Function
function [cycles,freqs] = createSimParams(c,f)
% CREATESIMPARAMS returns quantities needed by the Simulink model

    % Validate the input arguments
    arguments
        c {mustBeNumeric, mustBeVector} % Ensure array1 is numeric and a vector
        f {mustBeNumeric, mustBeVector} % Ensure array2 is numeric and a vector
        % buttonPushTime (3,1) double
    end

    % Check if both arrays have the same length
    if length(c) ~= length(f)
        error('Input arrays must be of the same length.');
    end

    cycles = c;
    freqs =  f;


end
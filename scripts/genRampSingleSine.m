function [signal_ALL,signal_ALL_time] = genRampSingleSine(freqVector,NumCycle_ramp,NumCycle_ss, TimeStep,CoolDownPeriod)
% User Input
freq = freqVector(:);
dt = TimeStep;

% Create Empty Cell
signal_stored = cell(size(freq));

for ifreq = 1:length(freq)

    %% Generate Trajectory
    NumCycle_total = 2*NumCycle_ramp + NumCycle_ss;
    OnePeriod = 1/freq(ifreq);
    TotalPeriod = OnePeriod*NumCycle_total;
    Nt = TotalPeriod/dt;
    t = ((1:Nt)-1)*dt;

    signal = sin(2*pi*freq(ifreq)*t);

    % Ramp signal
    WindLen = length(t);
    CosineFrac = 0.2;
    TukeyFilter = tukeywin(WindLen,CosineFrac);
    signal_ramped = TukeyFilter(:).*signal(:);

    signal_stored{ifreq} = signal_ramped;

end

% Concatenate signal
signal_ALL = [];
CoolDownSignal = zeros(1,CoolDownPeriod/dt);

for ifreq = 1:length(freq)
    signal_ALL = [signal_ALL(:); CoolDownSignal(:); signal_stored{ifreq}(:)];
end

% add cool down at the end
signal_ALL = [signal_ALL(:); CoolDownSignal(:)];

% create time variable
Nsignal_ALL = length(signal_ALL);
signal_ALL_time = ((1:Nsignal_ALL)-1)*dt;


end
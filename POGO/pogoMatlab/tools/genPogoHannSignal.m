function [ mOut ] = genPogoHannSignal( m, nCyc, freq, whichSig, whichFrame )
%genPogoHannSignal - add a Hann windowed signal to a model
%
%[ mOut ] = genPogoHannSignal( m, nCyc, freq[, whichSig[, whichFrame]] )
%
% m - input model
% mOut - output model
% nCyc, freq - number of cycles and frequency for signal

% whichSig - which signal number to apply it to (default 1)
% whichFrame - which frame it should be put in (default 1)
%
% m.dt and m.nt must be defined prior to running this.
%
%Written by P. Huthwaite, 2017
%Do not distribute.

if nargin < 5
    whichFrame = 1;
    if nargin < 4
        whichSig = 1;
    end
end



tSig = (0:m.nt-1)*m.dt;
sigLength = nCyc/freq;

wind = 0.5.*(1-cos(2*pi.*tSig./sigLength));
figure(10)
plot(wind)
wind(tSig < 0) = 0;
wind(tSig > sigLength) = 0;
figure(11)
plot(wind)


if isfield(m,'sigs') 
    error('Old version defining model.sigs is not permitted with this function. Use frames.')
end

m.frames{whichFrame}.sigs{whichSig}.sig = sin(2*pi*tSig*freq).*wind;

m.frames{whichFrame}.ntSig = m.nt;
m.frames{whichFrame}.dtSig = m.dt;

mOut = m;

end


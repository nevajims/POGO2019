function [ mOut ] = genPogoHannSignal_jim ( m, nCyc , freq , mode_shape , max_abs_sig  , whichSig, whichFrame )


%  function [ mOut ] = genPogoHannSignal( m, nCyc , freq, whichSig, whichFrame )
%  ---------------------------------------------------------------------------------------------------------------------------------------------------------------
%  ---------------------------------------------------------------------------------------------------------------------------------------------------------------
%  Jim Evans Mods-   modify the function so that it takes an additional parameter (complex number) and then uses it to apply a phase and amplitude to the signal 
%  mode_shape  ( ms_x,ms_y,ms_z )  -   
%  amp_ = abs(mode_shape);
%  phase_ = angle(mode_shape);
%  mode_shape is a complex number which give the phase and amplitude 
%  ---------------------------------------------------------------------------------------------------------------------------------------------------------------
%  ---------------------------------------------------------------------------------------------------------------------------------------------------------------


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
 wind(tSig < 0) = 0;
 wind(tSig > sigLength) = 0;


if isfield(m,'sigs') 
    error('Old version defining model.sigs is not permitted with this function. Use frames.')
end

% m.frames{whichFrame}.sigs{whichSig}.sig = sin(2*pi*tSig*freq).*wind;
amp_ = abs(mode_shape);
phase_ = angle(mode_shape);

m.frames{whichFrame}.sigs{whichSig}.sig =   (amp_* sin(2*pi*tSig*freq +phase_).*wind )/max_abs_sig;

% amp_z = abs(ms_z);
% phase_z = angle(ms_z);
% sig_z(:,node_index) = amp_z*sin(freq_range + phase_z).*wind;


m.frames{whichFrame}.ntSig = m.nt;
m.frames{whichFrame}.dtSig = m.dt;

mOut = m;



end


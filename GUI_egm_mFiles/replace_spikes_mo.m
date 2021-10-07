function [sig2,spk2] = replace_spikes_mo(signals,spikes,T,Th)
%% [sig2] = replace_spikes_mo(signals,spikes)
% This function detrend the electrograms and remove the spikes
%% IN %%
% - signals electrograms, as a matrix [LxN], L=length of the signals (ms); N = number of channels
% - spikes : temporal localization of spikes (ms)
% - T : % half width of a widow centered around the spike (interpolation in a window of width 2*T+1) [def. 15]
% - Th : Threshold to remove the spike [def. 10]
%% OUT %%
% - sig2 : detrended electrograms without the spikes
%%
% Michele Orini [12/2012]
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Default param
if nargin<4
    Th = 10;
end
if nargin<3
    T = 15;
end
if nargin<2
    error('This function requires signals and spikes as input')
end
%%
sig2 = nan(size(signals));
spk2=nan(length(spikes),size(signals,2));
% sigd = detrend(signals); % detrend
sigd = signals; 

for i = 1:size(signals,2)   % i=signals
    if sum(isnan(signals(:,i)))==0
        x=[];
        for j = 1:length(spikes) % j=heart beat
            W = spikes(j)+[-T:T]; % window centered around the spikes
            W(W<1)=[];W(W>length(sigd(:,i)))=[];
            tt = [W(abs(diff(sigd(W,i)))>Th*nanmedian(abs(diff(sigd(W,i)))))]; % samples to remove (those for which changes are too abrupt)
            if mean(diff(tt))>1 % to prevent no-continuous segments
                tt=tt(1):tt(end);
            end
            if ~isempty(tt)
                tt = [tt tt(end)+1];
                spk2(j,i) = tt(1);
            end
            x = [x(:);tt(:)];
        end
        sigd3 = sigd(:,i);sigd3(x)=[]; % remove samples
        t3=[1:size(sigd,1)];t3(x)=[];  % remove corresponding time instant
        sig2(:,i) = interp1(t3,sigd3,1:size(sigd,1),'linear','extrap'); % Interpolation
    end
end
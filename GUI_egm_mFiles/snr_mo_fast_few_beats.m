function [SNR,SNRdb,SNRdb_tv]=snr_mo_fast_few_beats(data,spikes,fs,Bsig,Bnoise,Nbeats);
% [SNR,SNRdb]=snr_mo_fast_few_beats(data,spikes,fs,Bsig,Bnoise);
% calclulate the signal-to-noise ratio in 20 beats (from DT to the end) and take the mean as
% estimate of SNR. In this way, signal contect is not over-estimated due to depolarization complexes and/or spikes, and harmonics of heart rate. 
% M Orini 04/2014

if nargin<6
    Nbeats = 1:length(spikes)-1;
end
if nargin<4
    Bsig = [2 45];
end
if nargin<5
    Bnoise = [45 100];
end

if sum(size(Bsig(:))==[2 1])+sum(size(Bnoise(:))==[2 1])~=4
    error('Bsig and Bnoise must be vectors containing 2 sclalars: lowest and uppest frequency (Hz) of signal and noise bands')
end
%
spikes = round(spikes/1000*fs);

% ic = randi(length(spikes)-1,[1,min(20,length(spikes))]);

% Find longer heart beats (to perform fft over more points)
Dsp = diff(spikes);
% if mean(Dsp > 1*fs)<0.3 % not to eliminate them in sinus R
% if ~length(Nbeats)==1
if length(Nbeats)>1
    if nanstd(Dsp)*fs/1000<10 % not to eliminate them in sinus R
        Dsp(Dsp > 1*fs) = nan;
    end
    [Dsps,isp] = sort(Dsp,'descend');
    if length(Nbeats(:))==1
        ic = isp(1:min(Nbeats,length(isp)-1));
    else
        ic = Nbeats;
    end
else
    ic = 1;
end
SNR = nan(size(data,2),length(ic));
SNRdb = nan(size(data,2),length(ic));
SNRdb_tv = nan(size(data,2),length(ic))';
% tic
for i = 1:length(ic)
    if length(ic)>1
    T = spikes(ic(i)) + round(120/1000*fs) : spikes(ic(i)+1)-20;
    T(T>size(data,1))=[];
    else
    T = spikes(ic(i)) + round(120/1000*fs) : min([size(data,1)-20,spikes(ic(i))+round(450/1000*fs)]);    
    end
    %     P = abs(fft(detrend(data(T,:)))).^2;
    if ~isnan(T)&~isempty(T)
        P = abs(fft((data(T,:)))).^2;
        f = [0:length(T)-1]/length(T)*fs;
        df = nanmedian(diff(f));
        Bs = [f>max(Bsig(1),df*1.2) & f<Bsig(2)];
        Bn = [f>Bnoise(1) & f<Bnoise(2)];
        SNR(:,i) = mean(P(Bs,:))./mean(P(Bn,:));
        SNRdb(:,i) = 10*log10(SNR(:,i));
        SNRdb_tv(ic(i),:) =  10*log10(SNR(:,i));
    end
end
% toc
% ic(ic>size(SNRdb,2))=[];

SNR = nanmean(SNR,2);
SNRdb = nanmean(SNRdb,2);

function [sig_corr,sig_corr_beats,sig_corr100,sig_corr100_beats,X] = sig_corr_mo(X,spikes,T,do_plot);
%% sig_corr = sig_corr_mo(X,spikes,do_plot);
% This function estimates the mean and the standard deviation of the intra-beat correlation, which is an idex of the regularity of the signal.
% IN :
% - X: X can be both a 2D and 3D matrix.
%         When it is a 2D matrix, X is interpreted as a matrix [TxNs] of signals, being T the length of the recordings and Ns the number of channels
%         When it is a 3D matrix, X is interpreted as a matrix [TxNbxNs] of signals, being T the length of the recordings, Nb the number of beats and Ns the number of channels
% - spikes: spikes is a vector whose elements are the time of stimulation (pacing) or activation
% - T : interval where the correlation is considered within each beat (samples)
% OUT:
% - sig_corr: sig_corr is a [Nsx2] matrix, containing mean and standard deviation of the correlation between each beats (waveform) and the average beat (waveform), for any channel/electrode
%% Michele Orini (01/2013)

%%
if nargin<4
    do_plot=0; % default
end
if nargin<3
    do_plot=0; % default
    T = [1:min(diff(spikes))];
end
if nargin<2&size(X,3)==1
    error('spikes, ie temporal localization of stimulation, is needed when X is a 2D matrix (array of EGM)')
end
%%

% if size(X,3)==1 % If input is a 2D array of EGM, reorganize it in 3D array [time,beats,channels]
%     signals_proc = X;
%     X = nan(min(diff(spikes)),length(spikes)-1,size(signals_proc,2));
%     for i=1:length(spikes)-1
%         X(1:min(diff(spikes)),i,:) = signals_proc(spikes(i):spikes(i)+min(diff(spikes))-1,:);
%     end
% end

%% Create matrix [time,heart beat,electrode] (just to save time)
if size(X,3)==1 % If input is a 2D array of EGM, reorganize it in 3D array [time,beats,channels]
    signals_proc = X;
    maxCL = T(end);
    L=min([max(diff(spikes)),maxCL]);
    X = nan(L,length(spikes),size(signals_proc,2));
    for i=1:length(spikes)-1
        H = spikes(i):spikes(i+1)-1;
        H(L:end)=[];H(H<1)=[];
            H(H>size(signals_proc,1))=[];

        X(1:length(H),i,:) = signals_proc(H,:);
    end
    H = spikes(i+1):size(signals_proc,1);
    H(L:end)=[];
    H(H>size(signals_proc,1))=[];
    X(1:length(H),i+1,:) = signals_proc(H,:);
end

sig_corr = nan(size(X,3),2);
sig_corr100 = nan(size(X,3),2);
sig_corr_beats= nan(size(X,2),2);
sig_corr100_beats= nan(size(X,2),2);

C = nan(size(X,2),size(X,3));
C100 = nan(size(X,2),size(X,3));
X2 = X;X2(isnan(X)) = 0;
T(T<1 | T>size(X2,1)) = [];
for i = 1:size(X,3)
    ci = corrcoef([nanmean(X2(T,:,i),2) X2(T,:,i)]); % correlation
    sig_corr(i,1) = nanmean(ci(2:end,1)); % mean of correlations between mean-waveform and the waveform of each beat
    sig_corr(i,2) = nanstd(ci(2:end,1)); % sd
    C(:,i) = ci(2:end,1);

    ci100 = corrcoef([nanmean(X2(100:end,:,i),2) X2(100:end,:,i)]); % correlation
    sig_corr100(i,1) = nanmean(ci100(2:end,1)); % mean of correlations between mean-waveform and the waveform of each beat
    sig_corr100(i,2) = nanstd(ci100(2:end,1)); % sd
    C100(:,i) = ci100(2:end,1);
end
sig_corr_beats(:,1) = nanmean(C,2); % mean of correlations between mean-waveform and the waveform of each beat
sig_corr_beats(:,2) = nanstd(C,1,2); % sd

sig_corr100_beats(:,1) = nanmean(C100,2); % mean of correlations between mean-waveform and the waveform of each beat
sig_corr100_beats(:,2) = nanstd(C100,1,2); % sd

if do_plot
    for i = 1:size(signals_proc,2)
        figure(10),
        subplot(211)
        hold off,plot(signals_proc(:,i));
        title(['#',num2str(i),' C=',num2str(Scorr_m(i),2),'\pm',num2str(Scorr_sd(i),2)])
        subplot(212)
        hold off,plot(X(:,:,i));
        pause
    end
end
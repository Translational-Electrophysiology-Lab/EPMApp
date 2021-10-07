function [sig_corr,sig_corr_m] = sig_corr_mo_all_beats(X,spikes,T);
%% [sig_corr,sig_corr_m] = sig_corr_mo_all_beats(X,spikes,T);
% = IN
% X : signals
% spikes : spikes. (samples)
% T : interval where to measure correlation (within one beat). (sample)
% = OUT
% sig_corr: Pearson's correlation between all beats [N*(N-1)/2,M] (N: number of beats, M: number of channels)
% sig_corr_m: Pearson's correlation between the average beat and the rest of beats [N,M] (N: number of beats, M: number of channels)
%%


if nargin<4
    do_plot=0; % default
end
if nargin<3
    do_plot=0; % default
    if length(spikes)>1
        T = [1:min(diff(spikes))];
    else
        T = [1:size(X,1)];
    end
end
    if nargin<2&size(X,3)==1
        error('spikes, ie temporal localization of stimulation, is needed when X is a 2D matrix (array of EGM)')
    end
    %%
    
    %% Create matrix [time,heart beat,electrode] (just to save time)
    if size(X,3)==1 % If input is a 2D array of EGM, reorganize it in 3D array [time,beats,channels]
        signals_proc = X;
        
        maxCL = T(end);
        L=min([max(diff(spikes)),maxCL]);
        X = nan(L,length(spikes),size(signals_proc,2));
        
        if length(spikes)>1
            for i= 1 : length(spikes)-1
                H = spikes(i) : spikes(i+1)-1;
                H(L:end)=[];H(H>size(signals_proc,1))=[];
                if H(1)<1
                    
                    X(-min(H)+H(H>0),i,:) =  signals_proc(H(H>0),:);
                else
                    H(H<1)=[];
                    X(1:length(H),i,:) = signals_proc(H,:);
                end
            end
            H = spikes(i+1):size(signals_proc,1);
            H(L:end)=[];
            X(1:length(H),i+1,:) = signals_proc(H,:);
        else
            i = 0;
            H = spikes(i+1):size(signals_proc,1);
            H(L:end)=[];H(H<1)=[];
            X(1:length(H),i+1,:) = signals_proc(H,:);
            
        end
    end
    
    sig_corr = nan(size(X,3),length(spikes)*(length(spikes)-1)/2);
    sig_corr_m = nan(size(X,3),1);
    X2 = X;X2(isnan(X)) = 0;
    T(T<1 | T>size(X2,1)) = [];
    for i = 1:size(X,3)
        ci = corrcoef([X2(T,:,i)]); % correlation between all beats
        a = triu(ci);%a(a==0|a==1)=[];
        a(isnan(a))=-10;
        a = a+tril(nan(size(a)));
        a = a(:);
        a(isnan(a))=[];
        a(a==-10)=nan;
        sig_corr(i,:) = a(:);
        
        ci2 = corrcoef([squeeze(nanmean(X2(T,:,i),2)) X2(T,:,i)]); % correlation between all beats
        sig_corr_m(i) = nanmean(ci2(1,2:end));
    end

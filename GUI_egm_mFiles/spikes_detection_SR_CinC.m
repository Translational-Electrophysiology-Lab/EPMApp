function [spikes,fwhm] = spikes_detection_SR_CinC(signals,fs,FWHM_thr);

if nargin<3
    FWHM_thr = inf; % ms
end

% sensthr_tot = [0.8 0.7 0.6 0.5 0.4 0.3 : -.025 : 0.05];
sensthr_tot = [0.8 0.6 0.5 0.4 0.3 0.2];

[a,b] = butter(3,[1 30]/(fs/2));
xflb = filtfilt(a,b,signals);
xf = medfilt1(diff(xflb).^2,round(50/1000*fs));
xf = xf./movmean(xf,round(10000/1000*fs)); % 09/04/2019
spikes_tot = cell(1,length(sensthr_tot));
W_tot = cell(1,length(sensthr_tot));
% rr_tot = cell(1,length(sensthr_tot));
for ithr = 1:length(sensthr_tot)
    clearvars -except FWHM_thr xflb fs signals ParamSig spikes_tot ithr sensthr_tot xf W_tot signals rr_tot rr_tot Xdir dir_load isbj filename dir_save
    sensthr = sensthr_tot(ithr);
    sp_all = cell(1,size(signals,2));
    sp_all_D = cell(1,size(signals,2));
    W = cell(1,size(signals,2)); % Duration of QRS complex
    
    %     h = xf>repmat(max(xf)*sensthr,[size(xf,1) 1]);
    h = xf>repmat(prctile(xf,99)*sensthr,[size(xf,1) 1]);
    h = imclose(h,ones(round(50/1000*fs),size(h,2)));
    %     xd = diff(signals);
    xd = diff(xflb);
    xd2 = diff(xd);
    for k = 1:size(h,2)
        ih = find(diff(h(:,k))>.99);
        eh = find(diff(h(:,k))<-.99); % modify to find VT/VF 02/02/2015
        
        if ~isempty(eh)&~isempty(ih)
            
            while eh(1)<ih(1)
                eh(1)=[];
                if isempty(eh)
                    break
                end
            end
            Lm = min(length(ih),length(eh));
            ih = ih(1:Lm); eh = eh(1:Lm);
            clear Lm
            % -
            sp_all{k} = nan(1,length(ih));
            W{k} = nan(1,length(ih));
            for j = 1:length(ih)
                H = (ih(j)-round(20/1000*fs)):(eh(j)+round(20/1000*fs)); % modify to find VT/VF 02/02/2015
                H(H<1 | H>size(xd2,1))=[];
                %                 ii= find([(xd(H(1:length(H)-1),k).*xd(H(2:length(H)),k))<=0 & xd2(H(1:end-1),k)<0])+H(1)+1;
                ii= find([(xd(H(1:length(H)-1),k).*xd(H(2:length(H)),k))<=0])+H(1)+1;
                x2 = detrend(signals(H,k),'constant');
                if length(ii)>1
                    %                 [~,id] = min(xd2(ii,k)); % min of second derivative: it
                    %                 works for clean and "spiky" QRS
                    
                    [~,id] = max(abs(x2(ii-H(1)))); %
                    ii = ii(id);
                end
                
                H2 = (ih(j)-round(60/1000*fs)):(eh(j)+round(60/1000*fs)); % modify to find VT/VF 02/02/2015
                H2(H2<1 | H2>size(signals,1))=[];
%                 FWHM = (eh(j)-ih(j))/fs*1000;
                FWHM = sum(abs(signals(H2,k))>[range(abs(signals(H2,k)))/2+min(abs(signals(H2,k)))])/fs*1000;
                if ~isempty(ii) && FWHM<=FWHM_thr
                    sp_all{k}(j) = ii-1;
                    W{k}(j) = FWHM;
                end
                
                %             figure,plot(H,xf(H,k)),hold on,plot(H,signals(H,k),'r'),plot(sp_all{k}(j),signals(sp_all{k}(j),k),'or')
            end
            iiko = sp_all{k}==0 | isnan(sp_all{k}) | [false diff(sp_all{k})<1];
            if ~isempty(iiko);
                sp_all{k}(iiko) = [];W{k}(iiko) = [];
            end
            sp_all_D{k} = diff(sp_all{k});
        end
        %     figure,plot(xf(:,k)),hold on,plot(signals(:,k),'r'),plot(sp_all{k}(:),signals(sp_all{k}(:),k),'or')
    end
    iiok = find(~cellfun(@isempty,sp_all_D));
    
    if ~isempty(iiok)
        [SD,iiSD] = min(cellfun(@nanstd,sp_all_D(iiok)));
        iiSD = iiok(iiSD);
        spikes = sp_all{iiSD};
        spikes = spikes/fs*1000;
        spikes_tot{ithr} = spikes;
        W_tot{ithr} = W{iiSD};
    else
        iiSD = [];
    end
end

% Eliminate configurations for which Q1 of RR < 200 ms (09/04/2019)
iiko = cellfun(@(x) prctile(diff(x),25),spikes_tot)<200;
spikes_tot(iiko) = [];
W_tot(iiko) = [];
spikes = [];fwhm = [];
if ~isempty(spikes_tot)
    %     [~,im]=min(cellfun(@(x) mad(diff(x,2),1)./nanmedian(diff(x,1)),spikes_tot));
    m = cellfun(@(x) mad(diff(x,2),0)./nanmedian(diff(x,1)),spikes_tot);
    ii = find(m < min(m)+eps*10 & m > min(m)-eps*10);
    if ~isempty(ii)
        spikes = spikes_tot{ii(end)};
        fwhm = W_tot{ii(end)};
    end
end

sp = round(spikes/1000*fs);
A = signals(sp,1);
if sum(abs(diff(sign(A))))>0
    w = mean(A(A>0))*mean(A>0)/abs( mean(A(A<0))*mean(A<0));
    if w<1
        for j = 1:length(sp)
            H = sp(j) + [-round((fwhm(j)+20)/1000*fs): round((fwhm(j)+20)/1000*fs)];
            H(H<1 | H>size(signals,1)) =[];
            [~,iM] = min(signals(H,1));
            sp(j) = H(iM);
        end
    else
        
        for j = 1:length(sp)
            H = sp(j) + [-round((fwhm(j)+20)/1000*fs): round((fwhm(j)+20)/1000*fs)];
            H(H<1 | H>size(signals,1)) =[];
            [~,iM] = max(signals(H,1));
            sp(j) = H(iM);
        end
        
    end
    
    spikes = sp/fs*1000;
end


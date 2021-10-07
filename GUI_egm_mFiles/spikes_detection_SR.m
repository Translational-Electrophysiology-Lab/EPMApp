function [spikes,iiSD] = spikes_detection_SR(signals,ParamSig);

sensthr_tot = [0.3 : -.025 : 0.05];
[a,b] = butter(3,[2 30]/(ParamSig.frequency/2));
xflb = filtfilt(a,b,signals);
xf = medfilt1(diff(xflb).^2,round(50/1000*ParamSig.frequency));
spikes_tot = cell(1,length(sensthr_tot));
rr_tot = cell(1,length(sensthr_tot));

for ithr = 1 : length(sensthr_tot)
    clearvars -except signals ParamSig spikes_tot ithr sensthr_tot xf signals rr_tot rr_tot Xdir dir_load isbj filename dir_save
    sensthr = sensthr_tot(ithr);
    sp_all = cell(1,size(signals,2));
    sp_all_D = cell(1,size(signals,2));
    
    %     h = xf>repmat(max(xf)*sensthr,[size(xf,1) 1]);
    h = xf>repmat(prctile(xf,98)*sensthr,[size(xf,1) 1]);
    h = imclose(h,ones(round(50/1000*ParamSig.frequency),size(h,2)));
    xd = diff(signals);
    xd2 = diff(xd);
    for k = 1:size(h,2)
        ih = find(diff(h(:,k))>.99);
        eh = find(diff(h(:,k))<-.99); % modify to find VT/VF 02/02/2015
        if ~isempty(eh)
            
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
            sp = nan(1,length(ih));
            for j = 1:length(ih)
                %             H = ih(j) + [-20:40];
                H = (ih(j)-round(20/1000*ParamSig.frequency)):(eh(j)+round(20/1000*ParamSig.frequency)); % modify to find VT/VF 02/02/2015
                H(H<1 | H>size(xd2,1))=[];
                ii= find([(xd(H(1:length(H)-1),k).*xd(H(2:length(H)),k))<=0 & xd2(H(1:end-1),k)<0])+H(1)+1;
                if length(ii)>1
                    %                 [~,id] = min(xd2(ii,k)); % min of second derivative: it
                    %                 works for clean and "spiky" QRS
                    x2 = detrend(signals(H,k));
                    [~,id] = max(abs(x2(ii-H(1)))); %
                    ii = ii(id);
                end
                if ~isempty(ii)
                    sp_all{k}(j) = ii-1;
                end
                
                %             figure,plot(H,xf(H,k)),hold on,plot(H,signals(H,k),'r'),plot(sp_all{k}(j),signals(sp_all{k}(j),k),'or')
            end
            sp_all{k}(sp_all{k}==0 | isnan(sp_all{k}) | [false diff(sp_all{k})<1]) = [];
            sp_all_D{k} = diff(sp_all{k});
        end
        %     figure,plot(xf(:,k)),hold on,plot(signals(:,k),'r'),plot(sp_all{k}(:),signals(sp_all{k}(:),k),'or')
    end
    iiok = find(~cellfun(@isempty,sp_all_D));
    
    if ~isempty(iiok)
        [SD,iiSD] = min(cellfun(@nanstd,sp_all_D(iiok)));
        iiSD = iiok(iiSD);
        spikes = sp_all{iiSD};
        spikes = spikes/ParamSig.frequency*1000;
        spikes_tot{ithr} = spikes;
        rr_tot{ithr} = diff(spikes,2);
    else
        iiSD = [];
    end
end

[~,im]=min(cellfun(@mad,rr_tot));
spikes = spikes_tot{im};
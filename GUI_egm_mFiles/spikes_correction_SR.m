function [spikes,ispikes] = spikes_correction_SR(signals,spikes,ParamSig,Niter,do_plot,ichan);

%% For Biobank
if nargin<6
    ichan = 1;
end

if nargin<5
    ichan = 1;
    
    Niter = 4;
end
if nargin<4
    ichan = 1;
    
    do_plot = 0;
    Niter = 4;
end

spikes_old = spikes;
spikes_samp_old = round(spikes_old/1000*ParamSig.frequency);
rr_old = diff(spikes_old);
spikes_tot = cell(1,Niter);
std_hrv_tot = nan(1,Niter);
%%
Thresh = 0.2;
Thresh_corr = 0.5;
MinBeat = 320; %ms
Temp_beat = 10; % number of beat considered to do the template
%%
maxQRS = round(120/1000*ParamSig.frequency);

for iter = 1:Niter
    %     if iter==1
    %         clearvars -except DataFromGui
    % %         do_plot=1;
    %         spikes = DataFromGui.spikes;
    %         signals = DataFromGui.signals;
    %         ParamSig = DataFromGui.ParamSig;
    %     else
    %         clearvars -except DataFromGui spikes signals ParamSig do_plot
    %     end
    clearvars -except spikes signals ParamSig do_plot iter spikes_tot *old Thresh* MinBeat Temp_beat maxQRS ichan std_hrv_tot
    if iter >1
        if isequal(spikes,spikes_tot{iter-1})
            
            [SDmin,ispikes] = min(std_hrv_tot);
            spikes = spikes_tot{ispikes};
            display(' * No further modifications')
            rr = diff(spikes);
            
            if do_plot
                figure,
                ax2(1)=subplot(211);
                plot(spikes_samp_old(2:end),rr_old,'.-k')
                hold on
                spikes_samp = round(spikes/1000*ParamSig.frequency);
                plot(spikes_samp(2:end),rr,'.--r')
                
                title(['SD(HRV) = ',num2str(nanstd(diff(spikes,2)),2),' ms; #iter=',num2str(iter-1)])
                ax2(2)=subplot(212);
                plot(signals(:,1)),hold on
                plot(spikes_samp,signals(spikes_samp,1),'xr','linewidth',2,'markersize',12)
                linkaxes(ax2,'x')
                xlabel('samples')
            end
            return
        end
    end
    
    spikes_samp = round(spikes/1000*ParamSig.frequency);
    rrm = nanmedian(diff(spikes));
    rrmd = diff(spikes,2);
    ii = find(abs(rrmd)>0.25*rrm)+1;
    rr = diff(spikes);
    Tmed = medfilt1(rr,20);
    %% Recognize Ectopics based on the average heart period and pattern (shorter-longer)
    %         Tm_ect = (rr(ii-1)+rr(ii))/2;
    %         for j = 1:length(Tm_ect)
    %             H = ii(j)+[-11:-2 3:12];
    %             H(H<1 | H>length(rr)) = [];
    %             Tm(j) = nanmean(rr(H));
    %             PatEct(j)= prod(sign(diff(spikes_samp(ii(j)-1:ii(j)+2),2)))<0;
    %         end
    %     ii_susp = ii( ~[abs(Tm_ect-Tm)./Tm<Thresh & PatEct]);
    %     % iiko = find(diff(ii_susp)==1);
    %     % ii_susp(iiko+1)=[];
    %     time_susp = spikes(ii_susp);
    %     Tm_ect_susp = Tm_ect(abs(Tm_ect-Tm)./Tm>Thresh);
    %     Tm_susp = Tm(abs(Tm_ect-Tm)./Tm>Thresh);
    %
    
    ii_susp = ii;
    time_susp = spikes(ii_susp);
    
    
    
    if do_plot
        
        figure
        plot(rr),
        hold on
        plot(ii,rr(ii),'or')
        plot(ii_susp,rr(ii_susp),'xk','markersize',12,'linewidth',2)
    end
    % figure,plot(signals(:,1)),hold on,
    % plot(spikes_samp,signals(spikes_samp,1),'ok'),
    % plot(spikes_samp([ii]),signals(spikes_samp([ii]),1),'og','markersize',8,'linewidth',2,'markerfacecolor','g'),
    % plot(spikes_samp([ii_susp]),signals(spikes_samp([ii_susp]),1),'xr','markersize',8,'linewidth',2),
    
    
    %% Create matrix of QRS [time,heart beat,electrode] (just to save time)
    DeltaT = round(30/1000*ParamSig.frequency);
    sp= spikes_samp -DeltaT; %samples
    Hend = DeltaT;
    L=min([max(diff(sp)),maxQRS]) + Hend;
    X = nan(L,length(sp),size(signals,2));
    for i=1:length(sp)-1
        H = sp(i):sp(i+1)-1+ Hend;
        H(L+1:end)=[];H(H<1)=[];H(H>length(signals))=[];
        X(1:length(H),i,:) = signals(H,:);
    end
    H = sp(i+1):size(signals,1);
    H(L:end)=[];
    X(1:length(H),i+1,:) = signals(H,:);
    clear sp
    %%
    Lb = zeros(size(ii_susp));
    sp2 = cell(1,length(ii_susp));
    
    for j = 1:length(ii_susp)
        ii_temp = ii_susp(j)+[-Temp_beat:-2 3:Temp_beat];
        ii_temp(ii_temp<1 | ii_temp>size(X,2)) = [];
        Xtemp = X(:,ii_temp,ichan);
        [cr,tau] = xcorr(Xtemp);
        [~,mm] = max(cr);
        Dm = tau(mm(1:size(Xtemp,2)));
        Xtemp_al = Xtemp;
        for j2 = 2:size(Xtemp,2)
            if Dm(j2)<0
                Xtemp_al(:,j2) = [Xtemp(-Dm(j2):end,j2);ones(-Dm(j2)-1,1)*Xtemp(end,j2)];
            elseif Dm(j2)>0
                Xtemp_al(:,j2) = [ones(Dm(j2)-1,1)*Xtemp(end,j2);Xtemp(1:end-Dm(j2)+1,j2)];
            end
        end
        clear tau cr Xtemp
        QRStemp = nanmean(Xtemp_al,2);
        clear H
        H = spikes_samp(ii_susp(j))-maxQRS:spikes_samp(ii_susp(j)+1)+maxQRS;
        H(H<1 | H>size(signals,1))= [];
        S = signals(H,ichan);
        clear H
        [cr,tau] = xcorr(S,QRStemp);
        iup = cr>(sum(QRStemp.^2)*Thresh_corr);
        iup = imclose(iup,ones(round(MinBeat/1000*ParamSig.frequency),1));
        iup = find(diff(iup)==1);
        
        if length(iup)>2
            iup2 = [];
            for j3 = 1:length(iup)
                H = iup(j3)+[0:round(100/1000*ParamSig.frequency)];
                H(H>length(cr))=[];
                [~,iup2(j3)] = max(cr(H));
                %                 H = iup(j3)+[0:round(100/1000*ParamSig.frequency)]-find(tau==0)+DeltaT;
                %                 H(H<1 | H>length(S))=[];
                %                 [~,iup2(j3)] = max(S(H));
            end
            sp2{j} = diff(iup(:)+iup2(:))-1;
            %             sp2{j} = diff(iup(:)+iup2(:));
            ie = abs(sp2{j}/ParamSig.frequency*1000-Tmed(ii_susp(j)))/(Tmed(ii_susp(j)))>Thresh;
            sp2{j}(ie)=[];
            iup(find(ie)+1)=[];
        end
        % To know which QRS I will eliminate later
        if length(iup)==1
            if tau(iup)>1.5*maxQRS
                iup=[];
            end
        end
        
        Lb(j) = length(iup);
        % %         figure,plot(QRStemp),hold on,plot(signals(spikes_samp(ii_susp(j))-maxQRS:spikes_samp(ii_susp(j)+1)+maxQRS ,ichan),'r')
        % %         figure,plot(tau,cr),hold on,plot([0 1000],[1 1]*(sum(QRStemp.^2)*.4))
        
        figure,
        tt = spikes_samp(ii_susp(j))-maxQRS:spikes_samp(ii_susp(j)+1)+maxQRS;
        tt(tt<1 | tt>size(signals,1))=[];
        plot(signals(:,ichan)),hold on
        plot(spikes_samp,signals(spikes_samp,ichan),'xk')
        plot(tt(1:length(QRStemp)),QRStemp,'k','linewidth',2),hold on,

        plot(tt,signals(tt,ichan),'-r')
        clear H
        plot(spikes_samp(ii_susp(j):ii_susp(j)+1),signals(spikes_samp(ii_susp(j):ii_susp(j)+1),ichan),'ok')
        if ~isempty(sp2(j))
            plot(spikes_samp(ii_susp(j))+cumsum(sp2{j}),signals(spikes_samp(ii_susp(j))+cumsum(sp2{j}),ichan),'xg','linewidth',3,'markersize',10)
        else
            plot(spikes_samp(ii_susp(j))+cumsum(sp2{j}),signals(spikes_samp(ii_susp(j))+cumsum(sp2{j}),ichan),'xg','linewidth',3,'markersize',10)
        end
        title(['(#B=',num2str(j),') Nb=',num2str(length(iup))])
        xlim([spikes_samp(ii_susp(j)) + [-1 1]*8*ParamSig.frequency]);
        close
    end
    clear j*
    %% OLD
    spikes_old = spikes;
    rr_old = rr;
    spikes_samp_old = spikes_samp;
    
    % if Lb==1 cancel the following bea, if Lb==0 cancel the first one
    spikes([ii_susp(Lb==1)+1 ii_susp(Lb==0)])=nan;
    %     spikes(ii_susp(Lb==1)+1)=nan;
    
    
    % add spikes
    % ii_add = find(Lb==3);
    % for j=1:length(ii_add)
    %     spikes = [spikes(1:ii_susp(ii_add(j))+(j-1)) (spikes_samp(ii_susp(ii_add(j)))+sp2{ii_add(j)}(1))/ParamSig.frequency*1000 spikes(ii_susp(ii_add(j))+1+(j-1) : end)];
    % end
    % ii_add = find(Lb>=3);
    % Nb = 0;
    % for j=1:length(ii_add)
    %     spikes = [spikes(1:ii_susp(ii_add(j))+Nb) (spikes_samp(ii_susp(ii_add(j)))+sp2{ii_add(j)}(1:end-1)')/ParamSig.frequency*1000 spikes(ii_susp(ii_add(j))+1+Nb : end)];
    %     Nb = Nb+length(sp2{ii_add(j)})-1;
    % end
    
    ii_add = find(Lb>=3);
    Nb = 0;
    spikes_add = [];
    for j=1:length(ii_add)
        spikes_add = [spikes_add (spikes_samp(ii_susp(ii_add(j))) + cumsum(sp2{ii_add(j)}(1:end-1))')/ParamSig.frequency*1000];
    end
    spikes = sort([spikes spikes_add]);
    spikes(isnan(spikes))=[];
    spikes(diff(spikes)==0)=[];
    spikes_samp = round(spikes/1000*ParamSig.frequency);
    rr = diff(spikes);
    
    if do_plot
        figure,
        ax(1)=subplot(211);
        plot(spikes_samp_old(2:end),rr_old/1000*ParamSig.frequency,'.-');
        hold on
        plot(spikes_samp(2:end),rr/1000*ParamSig.frequency,'o-r');
        
        ax(2)=subplot(212);
        plot(signals(:,1)),hold on
        plot(spikes_samp_old,signals(spikes_samp_old,1),'o')
        %     plot(spikes_samp([ii]),signals(spikes_samp([ii]),1),'og','markersize',5,'linewidth',5,'markerfacecolor','g'),
        plot(spikes_samp,signals(spikes_samp,1),'xr','linewidth',2,'markersize',12)
        plot(time_susp/1000*ParamSig.frequency,signals(round(time_susp/1000*ParamSig.frequency),1),'o','markersize',8,'linewidth',2,'markerfacecolor',[1 1 0])
        
        linkaxes(ax,'x')
        xlabel('samples')
    end
    spikes_tot{iter} = spikes;
    std_hrv_tot(iter) = std(diff(spikes,2));
    
    display(['Iter = ',num2str(iter), ': SD(HRV) = ',num2str(nanstd(diff(spikes,2)),4)])
    
end

[~,ispikes] = min(std_hrv_tot);
spikes = spikes_tot{ispikes};
rr = diff(spikes);

if do_plot
    figure,
    ax2(1)=subplot(211);
    plot(spikes_samp_old(2:end),rr_old,'.-k')
    hold on
    spikes_samp = round(spikes/1000*ParamSig.frequency);
    plot(spikes_samp(2:end),rr,'.--r')
    title(['SD(HRV) = ',num2str(nanstd(diff(spikes,2)),2),' ms; #iter=',num2str(iter)])
    
    ax2(2)=subplot(212);
    plot(signals(:,1)),hold on
    plot(spikes_samp,signals(spikes_samp,1),'xr','linewidth',2,'markersize',12)
    linkaxes(ax2,'x')
    xlabel('samples')
end

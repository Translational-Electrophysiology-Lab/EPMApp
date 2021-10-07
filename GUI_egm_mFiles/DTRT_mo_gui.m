function [Markers] = DTRT_mo_gui(signals_proc,spikes,ParamIn,DeltaT,types);
%% [Markers] = DTRT_mo_gui(signals_proc,spikes,ParamIn);

%% IN %%
% - signals_proc : electrograms (spikes already removed and band-pass filtered) as a matrix [LxN], L=length of the signals (ms); N = number of channels
% - spikes : temporal localization of spikes (ms)
% - ParamIn :
% - DeltaT : time before spike (in ms)
% - types : = 1 if EGM or ECG, = 2 if MAP
%% OUT %%
% - Markers: struct array with fields:
% dt: Activation Time (dV/dt)_min
% isot: Isopotential Time (not calculated for MAP)
% rt_up: Repolarization Time (dV/dt)_up [always with Wyatt method, only negative with alternative method]
% rt_down: Repolarization Time (dV/dt)_down [never with Wyatt method, only positive with alternative method]
% rt_Wyatt: Repolarization Time (dV/dt)_up (APD90 for MAPs)
% rt_Alternative: (dV/dt)_up for negative and (dV/dt)_down for positive T-waves
% ATw: Area under the T-wave: Area under V-V(isot) from isot to rt_up
% ATw2: Area under the T-wave: Area under V1-V1(1), with V1=detrend(V(isot:end))
% tdtM: Instant of V_max during activation
% tdtm: Instant of V_min during activation
% RatioI: [V(isot)-V(tdtm)]/[V(tdtM)-V(tdtm)]
% tTpeak : instant of T-Peak (Amplitude of MAP for MAP)
% ParamIn: parameters of analysis
%%
%% Example:
% ParamIn.DTmax =  150; % ms
% ParamIn.min_RT = 200; % ms
% ParamIn.max_RT = 480; % ms
% ParamIn.min_ARI = 150; % ms
% ParamIn.max_ARI = 380; % ms
% ParamIn.frequency = 2000; % Hz
% [Markers] = DTRT_mo_gui(signals_proc,spikes,ParamIn,10,1);

% michele orini 12/2012
% modified 02/2014
% modified 06/2014 (contraints on third derivative)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<4
    DeltaT = round(15/1000*ParamIn.frequency); % samples
    types = ones(size(signals_proc,2),1);
else
    if isempty(DeltaT)
        DeltaT = round(15/1000*ParamIn.frequency); % samples
    else
        DeltaT = round(DeltaT/1000*ParamIn.frequency);
    end
end

if nargin<5
    types = ones(size(signals_proc,2),1);
end
%% Param
if nargin<3|isempty(ParamIn)
    do_control=1;
    max_RT = round(nanmedian(diff(spikes))/3*2); %samples
    min_RT = round(200/1000*ParamIn.frequency); %samples
    max_ARI = round(nanmedian(diff(spikes))/3*2); %samples
    min_ARI = round(130/1000*ParamIn.frequency); %samples
    DTmax = round(180/1000*ParamIn.frequency); % samples
else
    if isfield(ParamIn,'do_control')
        do_control=ParamIn.do_control;
    else
        do_control=1;
    end
    max_ARI = round(ParamIn.max_ARI/1000*ParamIn.frequency); %samples
    min_ARI = round(ParamIn.min_ARI/1000*ParamIn.frequency); %samples
    max_RT = round(ParamIn.max_RT/1000*ParamIn.frequency); %samples
    min_RT = round(ParamIn.min_RT/1000*ParamIn.frequency); %samples
    DTmax = round(ParamIn.DTmax/1000*ParamIn.frequency); %samples
    
end
plotECM = 0;
%% this is just to consider that in few cases activation can occur before the spike (due to filtering)
% Note that this does not introduce any delay, because the temporal position of each marker does not change
% spikes= round( (spikes(:)-DeltaT)/1000*ParamIn.frequency); %samples
spikes= round( spikes(:)/1000*ParamIn.frequency -DeltaT); %samples
spikes(spikes<1)=1; % 23/01/2017

max_RT = max_RT + DeltaT;
min_RT = min_RT + DeltaT;
max_ARI = max_ARI + DeltaT;
min_ARI = min_ARI + DeltaT;
DTmax = DTmax + DeltaT;

Hend = 3*DeltaT;
%% Create matrix [time,heart beat,electrode] (just to save time)
maxCL = round(600/1000*ParamIn.frequency + DTmax);
L=min([max(diff(spikes)),maxCL]) + Hend;
X = nan(L,length(spikes),size(signals_proc,2));
for i=1:length(spikes)-1
    H = spikes(i):spikes(i+1)-1+ Hend;
    H(L:end)=[];H(H<1)=[];H(H>size(signals_proc,1))=[];
    X(1:length(H),i,:) = signals_proc(H,:);
end
H = spikes(end):size(signals_proc,1);
H(L:end)=[];
X(1:length(H),length(spikes),:) = signals_proc(H,:);

Xd = diff(X); Xd2 = diff(Xd); Xd3 = diff(Xd2);% Derivatives
% Initialization
tm = nan(size(X,2),size(X,3)); % dep time inside the heart beat
tiso = nan(size(X,2),size(X,3)); % iso time inside the heart beat
tMmin = nan(size(X,2),size(X,3));
tMmax= nan(size(X,2),size(X,3));
XdtM = nan(size(X,2),size(X,3));
Xdtm = nan(size(X,2),size(X,3));
tdtM = nan(size(X,2),size(X,3));
tdtm = nan(size(X,2),size(X,3));
Xiso = nan(size(X,2),size(X,3));
ATw= nan(size(X,2),size(X,3));
ATw2= nan(size(X,2),size(X,3));
iiTwpos = ones(1,size(X,3)); % if not able to identify I assume T-wave is positive
tTpeak = nan(size(X,2),size(X,3));
tTend_Tneg = nan(size(X,2),size(X,3));
tTend_Tpos = nan(size(X,2),size(X,3));
RatioI = nan(size(X,2),size(X,3));
% if sum(types==2)>0
%     tAmpMAP = nan(size(X,2),size(X,3));
%     AmpMAP = nan(size(X,2),size(X,3));
% end

%

Hwb = waitbar(0,'Calculating Markers ...');
%% Localize Markers
for i=1:size(X,3) % i=electrode
    waitbar(i/size(X,3),Hwb);
    
    if sum(isnan(signals_proc(:,i)))==0
        Ldt = DTmax+round(80/1000*ParamIn.frequency);
        ii = Xd2([1:Ldt],:,i).*Xd2([2:Ldt+1],:,i)<=0 & Xd([1:Ldt],:,i)<0 & Xd3([1:Ldt],:,i)>0; % zeros of II derivative & I derivative negative & III derivative > 0
        iim = Xd([1:Ldt],:,i).*Xd([2:Ldt+1],:,i)<0; % zeros of I derivative & I derivative negative
        iiMAP = Xd2([1:Ldt],:,i).*Xd2([2:Ldt+1],:,i)<0 & Xd([2:Ldt+1],:,i)>0 & Xd3([2:Ldt+1],:,i)<=0; % zeros of II derivative & I derivative positive & III derivative < 0(dV/dT max)
        if sum(ii(:))>0
            for j = 1:size(ii,2) % j=heart beat
                %% Depolarization time
                if types(i)==1
                    tt = find(ii(:,j));
                    if ~isempty(tt)
                        [~,kk]=min(Xd(tt,j,i));
                        tm(j,i) = tt(kk);
                        ttM = find(iim(:,j));
                        it1 = find(ttM<tt(kk) & Xd2(ttM,j,i)<0);
                        it2 = find(ttM>tt(kk) & ttM<tt(kk)+100 & Xd2(ttM,j,i)>0);
                        if ~isempty(it1)
                            it1=it1(end);
                            tdtM(j,i) = ttM(it1)+1;
                            XdtM(j,i) = X(ttM(it1),j,i);
                        else
                            tdtM(j,i) = 1;
                            XdtM(j,i) = X(1,j,i);
                        end
                        if ~isempty(it2)
                            it2=it2(1);
                            tdtm(j,i) = ttM(it2)+1;
                            Xdtm(j,i) = X(ttM(it2),j,i);
                        end
                    end
                elseif types(i)==2 % only for MAP
                    tt = find(iiMAP(:,j));
                    if ~isempty(tt)
                        [~,kk]=max(Xd(tt,j,i));
                        tm(j,i) = tt(kk);
                        
                        ttM = find(iim(:,j));%& Xd2([2:Ldt+1],j,i)<0);
                        it1 = find(ttM>tt(kk) & Xd2(ttM,j,i)<0);
                        if length(it1)>1
                            it1=it1(2);
                            tTpeak(j,i) = ttM(it1);
                        elseif length(it1)==1
                            it1=it1(1);
                            tTpeak(j,i) = ttM(it1);
                        end
                        
                    end
                else
                    error('types should be either =1 for EGM and ECG or =2 for MAP')
                end
            end
            clear ii iim Ldt
            %% iso time
            if types(i)==1 % only if signal is a ECG or EGM
                DTmaxiso = DTmax + round(100/1000*ParamIn.frequency);
                iidvdtmax = Xd2([1:DTmaxiso],:,i)<=0 & Xd([1:DTmaxiso],:,i)>=0; % II derivative negative & I derivative positive
                clear DTmaxiso
                if sum(iidvdtmax(:))>0
                    for j = 1:size(iidvdtmax,2) % j=heart beat
                        ttdvmax = find(iidvdtmax(:,j))+1;
                        if ~isempty(ttdvmax)
                            ttdvmax(ttdvmax<tm(j,i)|ttdvmax<tdtm(j,i)+15|ttdvmax < round(80/1000*ParamIn.frequency))=[];
                            if ~isnan(tm(j,i))
                                ttdvmax( [X(ttdvmax,j,i) > max(X(tm(j,i):end,j,i)) ] )=[];
                            else
                                ttdvmax = [];
                            end
                            if ~isempty(ttdvmax)
                                H = ttdvmax;
                                % To eliminate problems due to noise in depolarization phase
                                %                             ih = find(diff(H)>1);
                                %                             if ~isempty(ih) % eliminate small segments due to noise
                                %                                if ih<=15|ih/length(H)<.25
                                %                                    H(1:ih)=[];
                                %                                elseif length(H)-ih<=15
                                %                                    H(ih+1:end)=[];
                                %                                end
                                %                             end
                                if ~isempty(H)
                                    [~,tt]=min(Xd(H,j,i));
                                    tiso(j,i) = H(tt)-1;
                                    Xiso(j,i) = X(tiso(j,i),j,i);
                                end
                            end
                        end
                    end
                    % correction for XdtM: if tdtM==1 & XdtM<Xiso it means that the true XdtM occured before the spike
                    XdtM([tdtM(:,i)==1&XdtM(:,i)<Xiso(:,i)],i)=nan;
                end
                RatioI = (Xiso-Xdtm)./(XdtM-Xdtm);
            end
            %% Repolarization time (Wyatt+Alternative)
            W = min_RT:max_RT; % windows in which rep time can be estimated
            W(W>(size(Xd,1)-2))=[];
            ii = Xd2(W(1:end-1),:,i).*Xd2(W(2:end),:,i)<0 & Xd(W(2:end),:,i)>0 & Xd3(W(2:end),:,i)<=0; % zeros of II derivative & I derivative positive & III derivative < 0(dV/dT max)
            ii_alt = Xd2(W(1:end-1),:,i).*Xd2(W(2:end),:,i)<0 & Xd(W(2:end),:,i)<0 & Xd3(W(2:end),:,i)>=0;; % zeros of II derivative & I derivative positive & III derivative > 0 (dV/dT min)
            
            if types(i)==1
                for j=1:size(ii,2)
                    tt = find(ii(:,j))+W(1)-1;
                    % not too close to the spikes
                    %
                    if squeeze(nanmedian(sum(X(W,:,i)),2))>0 % Only for positive T-wave
                        if j<size(X,2)
                            tt(tt>spikes(j+1)-spikes(j)- round(10/1000*ParamIn.frequency))=[]; % Upstroke no too close to spike
                        end
                    end
                    tt(tt<=tiso(j,i))=[];
                    tt(tt>tm(j,i)+max_ARI|tt<tm(j,i)+min_ARI) = []; % control over ARI
                    if ~isempty(tt)
                        if squeeze(nanmedian(sum(X(W,:,i)),2))>0 % Only for positive T-wave
                            if ~isempty(max(find([Xd2(tt,j,i)-Xd2(tt-1,j,i)]>0)))&max(find([Xd2(tt,j,i)-Xd2(tt-1,j,i)]>0))<length(tt)& max(find([Xd2(tt,j,i)-Xd2(tt-1,j,i)]>0))>1
                                tt(1:max(find([Xd2(tt,j,i)-Xd2(tt-1,j,i)]>0)))=[]; % eliminate fiducial points which occurr at same time or before the first change from convex to concave (can be part of a slow depolarization)
                            end
                        end
                        [p1up,ii2] = max(Xd(tt,j,i));
                        tMmax(j,i)=tt(ii2);
                    end
                    % downslope
                    ttalt = find(ii_alt(:,j))+W(1)-1;
                    if j<size(ii,2)
                        ttalt(ttalt<tMmax(j,i) | spikes(j+1)-spikes(j)-ttalt < round(10/1000*ParamIn.frequency))=[];
                    else
                        ttalt(ttalt<tMmax(j,i))=[];
                    end
                    if ~isempty(ttalt)
                        [p1down,ii2] = min(Xd(ttalt,j,i));
                        tMmin(j,i)=ttalt(ii2);
                    end
                    %% T-end
                    % Measuring T-end in both positive and negative Twave, I will decide later
                    if ~isnan(tMmin(j,i))
                        r = roots([p1down  (X(tMmin(j,i),j,i)-p1down*tMmin(j,i))]);
                        if r<size(X,1) % no longer than cycle length
                            tTend_Tpos(j,i) = r; % time for which the line tangent to tMin is equal to zero (for positive T waves)
                        end
                    end
                    
                    if ~isnan(tMmax(j,i))
                        r = roots([p1up  (X(tMmax(j,i),j,i)-p1up*tMmax(j,i))]);
                        if r<size(X,1) % no longer than cycle length
                            tTend_Tneg(j,i) = r; % time for which the line tangent to tMax is equal to zero (for positive T waves)
                        end
                    end
                    %%
                    % Area Twave
                    tr = max([tMmax(j,i),tMmin(j,i),round(max_RT*.8)]);
                    if ~isnan(tiso(j,i)+tr); % mo 17/06/2015
                        Xbl = X(:,j,i)-X(tiso(j,i),j,i);
                        tr = min([tr,size(Xbl,1)]);
                        ATw(j,i) = sum(Xbl(tiso(j,i):tr))/nansum(abs(Xbl(tiso(j,i):tr)));
                        
                    end
                    % Area Twave2
                    if ~isnan(tiso(j,i));
                        A = X(tiso(j,i):end,j,i);A(isnan(A))=[];
                        Xbl = detrend(A);clear A
                        Xbl=Xbl-Xbl(1);
                        ATw2(j,i) = sum(Xbl)/nansum(abs(Xbl));
                    end
                    
                    
                end
            elseif types(i)==2 % MAP
                ATw2(:,i) = 10;
                ATw(:,i) = 10;
                for j = 1:size(ii,2)
                    if ~isnan(tTpeak(j,i))
                        H = find(X(W,j,i)< (.1*X(tTpeak(j,i),j,i) + prctile(X(:,j,i),10))); % prctile(X(:,j,i),10)) represents repolarization voltages
                        if ~isempty(H)
                            tMmax(j,i) = H(1)+W(1); % MAP90 as Wyatt, (dV/dt)min as alternative
                        end
                    end
                    % dv/dt-min
                    ttalt = find(ii_alt(:,j))+W(1)-1;
                    if j<size(ii,2)
                        ttalt(spikes(j+1)-spikes(j)-ttalt < round(10/1000*ParamIn.frequency))=[];
                    else
                        ttalt(ttalt<tMmax(j,i))=[];
                    end
                    if ~isempty(ttalt)
                        [p1down,ii2] = min(Xd(ttalt,j,i));
                        tMmin(j,i)=ttalt(ii2);
                    end
                    
                end
            else
                error('types should be either =1 for EGM and ECG or =2 for MAP')
                
            end
            
            
            %% polarity of Twave (per channel)
            %% Rules based on Ratio only works for paced UEG
            
%             if mean(isnan(ATw2(:,i)))<1&mean(isnan(RatioI(:,i)))<.5
%                 if abs(nanmean(ATw2(:,i)))>.80 % if enough positive/negative => rely on area only
%                     iiTwpos(i)=sign(nanmean(ATw2(:,i)))>0;
%                 else
%                     iiTwpos(i) = [[nanmedian(ATw2(:,i))>=0 | nanmedian(RatioI(:,i))>0.80]  & nanmedian(RatioI(:,i))>0.25];
%                 end
%             elseif mean(isnan(ATw(:,i)))<1&mean(isnan(RatioI(:,i)))<.5
%                 if abs(nanmean(ATw(:,i)))>.80
%                     iiTwpos(i)=sign(nanmean(ATw(:,i)))>0;
%                 else
%                     iiTwpos(i) = [[nanmedian(ATw(:,i))>=0 | nanmedian(RatioI(:,i))>0.80]  & nanmedian(RatioI(:,i))>0.25];
%                 end
%             elseif mean(isnan(ATw2(:,i)))<1&mean(isnan(RatioI(:,i)))>.5
%                 iiTwpos(i) = nanmedian(ATw2(:,i))>=0;
%             elseif mean(isnan(ATw(:,i)))<1&mean(isnan(RatioI(:,i)))>.5
%                 iiTwpos(i) = nanmedian(ATw(:,i))>=0;
%             end
 
            iiTwpos(i) = nanmedian(ATw(:,i))>=0;

            %% T-peak
            W = [min_RT:max_RT]; % windows in which rep time can be estimated
            W(W>(size(Xd,1)-1))=[];
            iiMax = Xd(W(1:end-1),:,i).*Xd(W(2:end),:,i)<0 & Xd2(W(2:end),:,i)<0; % zeros of I derivative & II derivative negative (V max)
            iiMin = Xd(W(1:end-1),:,i).*Xd(W(2:end),:,i)<0 & Xd2(W(2:end),:,i)>0; % zeros of I derivative & II derivative positive (V min)
            if types(i)==1
                for j=1:size(ii,2)
                    if iiTwpos(i)
                        tt = W(find(iiMax(:,j)))+1;
                        tt(tt<=tMmax(j,i)|tt<tiso(j,i)|tt>=tMmin(j,i))=[];
                        if ~isempty(tt)
                            [~,i2]=max(X(tt,j,i));
                            tTpeak(j,i)=tt(i2);
                        end
                        clear tt i2
                    else% negative Twave
                        tt = W(find(iiMin(:,j)))+1;
                        tt(tt>=tMmax(j,i)|tt<tiso(j,i))=[];
                        if ~isempty(tt)
                            [~,i2]=min(X(tt,j,i));
                            tTpeak(j,i)=tt(i2);
                            clear tt i2
                        end
                        
                    end
                end
            end
        end
        %% plot control
        %                 figure,plot(X(:,j,i)),hold on,
        %                 plot([0 600],[1 1]*prctile(X(1:tiso(j,i),j,i),75),'--k')
        %                 plot([0 600],[1 1]*prctile(X(1:tiso(j,i),j,i),25),'--k')
        %                 plot([0 600],[1 1]*prctile(X(1:tiso(j,i),j,i),66),'--g')
        %                 plot([0 600],[1 1]*prctile(X(1:tiso(j,i),j,i),33),'--g')
        %                 plot(tiso(j,i),X(tiso(j,i),j,i),'squarer'),
        %                 plot(tMmin(j,i),X(tMmin(j,i),j,i),'or'),
        %                 plot(tMmax(j,i),X(tMmax(j,i),j,i),'xr','markersize',10,'linewidth',2),
        %                 plot(tMalt(j,i),X(tMalt(j,i),j,i),'+g','markersize',10,'linewidth',2)
        
        %
        % %
        % % %
        
        if plotECM
            figure(10),
            ax(1)=subplot(311);
            hold off
            plot(signals_proc(:,i)),hold on
            dtplot=tm(:,i)+spikes(1:size(tm,1));dtplot(isnan(dtplot))=[];
            rtplot=tMmax(:,i)+spikes(1:size(tMmax,1));rtplot(isnan(rtplot))=[];
            rtplot2=tMmin(:,i)+spikes(1:size(tMmin,1));rtplot2(isnan(rtplot2))=[];
            isoplot=tiso(:,i)+spikes(1:size(tMmax,1));isoplot(isnan(isoplot))=[];
            tdtmplot = tdtm(:,i)+spikes(1:size(tMmax,1));tdtmplot(isnan(tdtmplot))=[];
            tdtMplot = tdtM(:,i)+spikes(1:size(tMmax,1));tdtMplot(isnan(tdtMplot))=[];
            plot(dtplot,signals_proc(dtplot,i),'or')
            plot(rtplot,signals_proc(rtplot,i),'+r','markersize',8,'linewidth',2)
            plot(rtplot2,signals_proc(rtplot2,i),'xg','markersize',8,'linewidth',2)
            plot(isoplot,signals_proc(isoplot,i),'squarer')
            plot(tdtmplot,signals_proc(tdtmplot,i),'^k')
            plot(tdtMplot,signals_proc(tdtMplot,i),'vk')
            for jj = 1:length(spikes)
                hold on,plot([1 1]*spikes(jj),get(gca,'ylim'),'--k','linewidth',2)
            end
            %             for jj = 1:length(spikes)
            %                 hold on,plot([1 1]*nanmedian(tm(:,jj))+spikes(jj)+min_ari,get(gca,'ylim'),'--k')
            %                 hold on,plot([1 1]*nanmedian(tm(:,jj))+spikes(jj)+max_ari,get(gca,'ylim'),'--k')
            %             end
            title(['electrode #',num2str(i)])
            %%
            ax(2)=subplot(312);
            hold off
            plot(tMmax(:,i)+spikes(1:size(tMmax,1)),tMmax(:,i),'.-r'),hold on
            plot(tMmin(:,i)+spikes(1:size(tMmin,1)),tMmin(:,i),'.-g'),hold on
            plot(tm(:,i)+spikes(1:size(tm,1)),tm(:,i),'.-k')
            title('Repolarization & depolarization time')
            legend('Rep(up)','Rep(down)','Dep')
            xlabel('Time [ms]'),ylabel('Time [ms]')
            %%
            ax(3)=subplot(313);
            hold off
            plot(tMmax(:,i)+spikes(1:size(tMmax,1)),tMmax(:,i)-tm(:,i),'.--r')
            hold on
            plot(tMmin(:,i)+spikes(1:size(tMmin,1)),tMmin(:,i)-tm(:,i),'.--g')
            legend('ARIup','ARIdown')
            title(['ARI  with RTup and RTdown']),xlabel('Time [ms]'),ylabel('Time [ms]')
            
            linkaxes(ax,'x')
            pause
            
            
        end
    end
end

dt = tm+repmat(spikes,[1 size(X,3)]) ; % unwrap
tdtM(XdtM<0|isnan(XdtM))=nan;
ic = mean(isnan(tdtM))>0.5;
tdtM(:,ic)=nan;clear ic
tRw = tdtM+repmat(spikes,[1 size(X,3)]) ; % unwrap
tSw = tdtm+repmat(spikes,[1 size(X,3)]) ; % unwrap
rt_up = tMmax+repmat(spikes,[1 size(X,3)]) ; % Unwrap
rt_down = tMmin+repmat(spikes,[1 size(X,3)]) ; % Unwrap
isot = tiso+repmat(spikes,[1 size(X,3)]) ; % unwrap
tTpeak = tTpeak+repmat(spikes,[1 size(X,3)]) ; % unwrap
tTend_Tpos = tTend_Tpos+repmat(spikes,[1 size(X,3)]) ; % unwrap
tTend_Tneg = tTend_Tneg+repmat(spikes,[1 size(X,3)]) ; % unwrap

if do_control
    %% control DT
    dtold = dt;
    tmold = tm;
    varDep = nanstd(diff(dt));
    % iidep = find(varDep>3*nanmean(varDep));%iidep(SNR(iidep)<20)=[];
    iidep = find(varDep>2.5*nanmedian(varDep));%iidep(SNR(iidep)<20)=[];
    for i = 1:length(iidep)
        ie = iidep(i);
        C = corrcoef([ squeeze(nanmean(X(1:min([2*max(dt(:,ie)-spikes),200]),:,ie),2)) X(1:min([2*max(dt(:,ie)-spikes),200]),:,ie)]);
        C = C(2:end,1)';
        dtdum = dt(:,ie)-spikes;
        C(dtdum<=10)=nan;dtdum(dtdum<=5)=nan; % only dt>10 ms
        dtm = round(nanmedian(dtdum(C>.9)));  % only dt>10 ms & C>0.9
        if ~isnan(dtm)
            WD = dtm+[-20:20];WD(WD<2)=[];
            ii = Xd2(WD,:,ie).*Xd2(WD-1,:,ie)<0 & Xd(WD,:,ie)<0;
            tm2 = nan(size(ii,2),1);
            for j = 1:size(ii,2)
                tt = find(ii(:,j));
                if ~isempty(tt)
                    [~,kk]=min(abs(tt-dtm));
                    tm2(j) = tt(kk);
                end
            end
            tm(:,ie) = tm2+(dtm-22);
        end
        %    figure,plot(tmold(:,ie),'.-'),hold on,plot(tm(:,ie),'.-r')
        %    legend('DTold','DTnew')
        %    title(num2str(ie))
        %    pause
        % --
    end
    dt = tm+repmat(spikes,[1 size(X,3)]);
    
    %% control RT
    rtold = rt_up;
    tMold = tMmax;
    varRep = nanstd(diff(rt_up));
    % iirep = find(varRep>3*nanmean(varRep));%iirep(SNR(iirep)<20)=[];
    iirep = find(varRep>2.5*nanmedian(varRep));%iirep(SNR(iirep)<20)=[];
    
    for i = 1:length(iirep)
        ie = iirep(i);
        %         t0 = round(nanmedian(tmold(:,ie))) + min_ari;
        t0 = min_RT;
        C = corrcoef([ squeeze(nanmean(X(t0:end,:,ie),2)) X(t0:end,:,ie)]);
        C = C(2:end,1);
        rtdum = rtold(:,ie)-spikes;
        rtm = round(nanmedian(rtdum(C>.90)));
        ii = Xd2(t0:end-1,:,ie).*Xd2(t0+1:end,:,ie)<0 & Xd(t0+2:end,:,ie)>0;
        tM2 = nan(size(ii,2),1);
        for j = 1:size(ii,2)
            tt = find(ii(:,j))+t0-1;
            if ~isempty(tt)
                [~,kk]=min(abs(tt-rtm));
                tM2(j) = tt(kk);
            end
        end
        tMmax(:,ie) = tM2;
        
    end
    rt_up = tMmax+repmat(spikes,[1 size(X,3)]);
end

%%
leg = {'dt [ms]: Activation Time (dV/dt)_min',...,
    'isot [ms]: Isopotential Time',...,
    'rt_up [ms]: Repolarization Time (dV/dt)_up [always with Wyatt method, only negative with alternative method]',...,
    'rt_down [ms]: Repolarization Time (dV/dt)_down [never with Wyatt method, only positive with alternative method]',...,
    'rt_Wyatt [ms]: Repolarization Time (dV/dt)_up',...,
    'rt_Alternative [ms]: (dV/dt)_up for negative and (dV/dt)_down for positive T-waves',...,
    'tTpeak [ms]: instant of T-Peak',...,
    'ATw: Area under the T-wave: Area under V-V(isot) from isot to rt_up',...,
    'ATw2: Area under the T-wave: Area under V1-V1(1), with V1=detrend(V(isot:end))',...,
    'tdtM [ms]: Instant of V_max during activation',...,
    'tdtm [ms]: Instant of V_min during activation',...,
    'tRw [ms]: R-wave (as tdtM)',...,
    'tSw [ms]: S-wave (as tdtm)',...,
    'RatioI: [V(isot)-V(tdtm)]/[V(tdtM)-V(tdtm)]',...,
    'ParamIn: parameters of analysis'};


iiTwpos = logical(iiTwpos); % true if T-wave is positive, false otherwise
% from samples to milli-seconds
Markers.tTpeak = tTpeak/ParamIn.frequency*1000;

Markers.tTend = nan(size(dt));
Markers.tTend(:,iiTwpos) = tTend_Tpos(:,iiTwpos)/ParamIn.frequency*1000;
Markers.tTend(:,~iiTwpos) = tTend_Tneg(:,~iiTwpos)/ParamIn.frequency*1000;
% Markers.tTend(end,[Markers.tTend(end,:)>size(signals_proc,1)]) = nan;
Markers.tTend(Markers.tTend/1000*ParamIn.frequency>size(signals_proc,1)) = nan;


Markers.rt_Wyatt = rt_up/ParamIn.frequency*1000;
Markers.rt_Alternative(:,~iiTwpos) =rt_up(:,~iiTwpos)/ParamIn.frequency*1000;
Markers.rt_Alternative(:,iiTwpos) = rt_down(:,iiTwpos)/ParamIn.frequency*1000;

Markers.tSw = tSw/ParamIn.frequency*1000;
Markers.tRw = tRw/ParamIn.frequency*1000;
Markers.iiTwpos = iiTwpos;
Markers.dt = dt/ParamIn.frequency*1000;
Markers.rt_up = rt_up/ParamIn.frequency*1000;
Markers.isot=isot/ParamIn.frequency*1000;
Markers.rt_down=rt_down/ParamIn.frequency*1000;
Markers.ATw=ATw;
Markers.ATw2=ATw2;
Markers.tdtM = tdtM/ParamIn.frequency*1000;
Markers.tdtm = tdtm/ParamIn.frequency*1000;
Markers.RatioI=RatioI;
Markers.ParamIn = ParamIn;
Markers.m_function = {'DTRT_mo_gui.m'};
Markers.Legend = leg;

close(Hwb)



% Ex DTRT_mo_gui_bio_bank_v3.m
% Same as DTRT_mo_gui_v2.m but with only 1 for circle

function [Markers,Hwb] = DTRT_mo_gui_v4(signals_proc,spikes,ParamIn,Info_Correct);
%% [Markers] = DTRT_mo_gui_bio_bank(signals_proc,spikes,ParamIn,DeltaT,Tw_polarity);
% I have deleted ARI_max and ARI_min and made RT_min and RT)max time-varying

%% IN %%
% - signals_proc : electrograms (spikes already removed and band-pass filtered) as a matrix [LxN], L=length of the signals (ms); N = number of channels
% - spikes : temporal localization of spikes (ms)
% - ParamIn :
% - DeltaT : time before spike (in ms)
%% OUT %%
% - Markers: struct array with fields:
% see legend
%%
%% Example:
% ParamIn.DTmax =  150; % ms
% ParamIn.min_RT = 200; % ms
% ParamIn.max_RT = 480; % ms
% ParamIn.frequency = 2000; % Hz
% [Markers] = DTRT_mo_gui(signals_proc,spikes,ParamIn,10,1);

% michele orini 12/2012
% modified 02/2014
% modified 06/2014 (contraints on third derivative)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Param
if nargin<4
    Info_Correct = [];
end
if nargin<3|isempty(ParamIn);
    Info_Correct = [];
    answ = inputdlg({'Sampling frequency (Hz)','AT max (ms)','RT min (ms)','RT max (ms)','DeltaT (ms)'},'Inputs',1,{'1000','150','200','400','15'});
    if ~isempty(answ);
        ParamIn.frequency = str2double(answ{1});
    else
        return
    end
    max_RT = round(str2double(answ{4})/1000*ParamIn.frequency); %samples
    min_RT = round(str2double(answ{3})/1000*ParamIn.frequency); %samples
    DTmax = round(str2double(answ{2})/1000*ParamIn.frequency); % samples
    DeltaT = round(str2double(answ{5})/1000*ParamIn.frequency); %samples
else
    max_RT = round(ParamIn.max_RT/1000*ParamIn.frequency); %samples
    min_RT = round(ParamIn.min_RT/1000*ParamIn.frequency); %samples
    DTmax = round(ParamIn.DTmax/1000*ParamIn.frequency); %samples
    DeltaT = round(ParamIn.DeltaT/1000*ParamIn.frequency); % samples
end


%%
% time-varying Repol. boundaries
if numel(max_RT)==1
    max_RT = ones(1,length(spikes))*max_RT;
end
if numel(min_RT)==1
    min_RT = ones(1,length(spikes))*min_RT;
end

%%
tPi_max = 230/1000*ParamIn.frequency; %samples
tPi_min = 80/1000*ParamIn.frequency; %samples
QRSw_max = 200/1000*ParamIn.frequency; %samples

%% this is just to consider that in few cases activation can occur before the spike (due to filtering)
% Note that this does not introduce any delay, because the temporal position of each marker does not change
% spikes= round( (spikes(:)-DeltaT)/1000*ParamIn.frequency); %samples
spikes= round( spikes(:)/1000*ParamIn.frequency -DeltaT); %samples


max_RT = max_RT + DeltaT;
min_RT = min_RT + DeltaT;
DTmax = DTmax + DeltaT;

if ~isempty(Info_Correct)
    Nchannel = Info_Correct.ichan;
    signals_proc = signals_proc(:,Nchannel);
else
    Nchannel = 1:size(signals_proc,2);
end

Hend = 3*DeltaT;
%% Create matrix [time,heart beat,electrode] (just to save time)
maxCL = round(600/1000*ParamIn.frequency + DTmax);
% L = min([max(diff(spikes)),maxCL]) + Hend;
L = maxCL;

X = nan(L,length(spikes),size(signals_proc,2));
if length(spikes)>1
    for i= 1 : length(spikes)-1
        H = spikes(i) : spikes(i+1)-1+Hend;
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
    i = 1;
    H = spikes(i):size(signals_proc,1);
    H(L:end)=[];
    if H(1)<1
        X(-min(H)+H(H>0),i,:) =  signals_proc(H(H>0),:);
    else
        H(H<1)=[];
        X(1:length(H),i,:) = signals_proc(H,:);
    end
end

%%
if ~isempty(Info_Correct)
    Nbeats = Info_Correct.Nbeats;
    max_RT = max_RT(Info_Correct.Nbeats);
    min_RT = min_RT(Info_Correct.Nbeats);
else
    Nbeats = 1:size(X,2);
end
% -
% if ~isempty(Info_Correct)
%     Nchannel = Info_Correct.ichan;
%     signals_proc = signals_proc(:,Nchannel);
% else
%     Nchannel = size(X,3);
% end

%
X = X(:,Nbeats,:);
Xd = diff(X); Xd2 = diff(Xd); Xd3 = diff(Xd2);% Derivatives


%% Initialization
if isempty(Info_Correct);
    tm = nan(size(X,2),size(X,3)); % dep time inside the heart beat
    tiso = nan(size(X,2),size(X,3)); % iso time inside the heart beat
    tMmin = nan(size(X,2),size(X,3));
    tMmax= nan(size(X,2),size(X,3));
    XdtM = nan(size(X,2),size(X,3));
    Xdtm = nan(size(X,2),size(X,3));
    tdtM = nan(size(X,2),size(X,3)); % R-wave (defined as the prominent wave, independently of polarity)
    tdtm = nan(size(X,2),size(X,3));
    Xiso = nan(size(X,2),size(X,3));
    ATw= nan(size(X,2),size(X,3));
    Tpeak_amp = nan(size(X,2),size(X,3));
    ATw2= nan(size(X,2),size(X,3));
    iiTwpos = ones(1,size(X,3)); % if not able to identify I assume T-wave is positive
    tTpeak = nan(size(X,2),size(X,3));
    tTend_Tneg = nan(size(X,2),size(X,3));
    tTend_Tpos = nan(size(X,2),size(X,3));
    RatioI = nan(size(X,2),size(X,3));
    tPw = nan(size(X,2),size(X,3));
    QRSw_fw90 = nan(size(X,2),size(X,3)); % QRSwith robust (as full width at 90& of duration)
    tQw = nan(size(X,2),size(X,3)); % Instant of Q wave (defined as the deflection preceeding the prominent one)
    tSw = nan(size(X,2),size(X,3)); % Instant of S wave (defined as the deflection following the prominent one)
    tQRSon = nan(size(X,2),size(X,3)); % Flex point before Q wave
    tQRSoff = nan(size(X,2),size(X,3)); % Flex point after S wave
    tTend_defl = nan(size(X,2),size(X,3)); % Tend as the deflection point following the Tpeak
    
    Rw_amp = nan(size(X,2),size(X,3));
    Sw_amp = nan(size(X,2),size(X,3));
    Qw_amp = nan(size(X,2),size(X,3));
    QRS_amp = nan(size(X,2),size(X,3)); % peak to peak (range(QRS))
    QRS_area = nan(size(X,2),size(X,3)); % Area from Qw to Sw
    
else
    
    if isfield(Info_Correct.MarkersC,'dt');
        tm = round( (Info_Correct.MarkersC.dt(Nbeats,Nchannel) - spikes(Nbeats)*ones(1,length(Nchannel))/ParamIn.frequency*1000)/1000*ParamIn.frequency);
    else
        tm = nan(size(X,2),size(X,3));
    end
    
    if isfield(Info_Correct.MarkersC,'tiso');
        tiso = round( (Info_Correct.MarkersC.isot(Nbeats,Nchannel) - spikes(Nbeats)*ones(1,length(Nchannel))/ParamIn.frequency*1000)/1000*ParamIn.frequency); % iso time inside the heart beat
    else
        tiso = nan(size(X,2),size(X,3));
    end
    
    if isfield(Info_Correct.MarkersC,'rt_down');
        tMmin = round( (Info_Correct.MarkersC.rt_down(Nbeats,Nchannel) - spikes(Nbeats)*ones(1,length(Nchannel))/ParamIn.frequency*1000)/1000*ParamIn.frequency);
    else
        tMmin = nan(size(X,2),size(X,3));
    end
    
    if isfield(Info_Correct.MarkersC,'rt_up');
        tMmax= round( (Info_Correct.MarkersC.rt_up(Nbeats,Nchannel) - spikes(Nbeats)*ones(1,length(Nchannel))/ParamIn.frequency*1000)/1000*ParamIn.frequency);
    else
        tMmax = nan(size(X,2),size(X,3));
    end
    
    if isfield(Info_Correct.MarkersC,'XdtM');
        XdtM= Info_Correct.MarkersC.Rw_amp(Nbeats,Nchannel);
    else
        XdtM = nan(size(X,2),size(X,3));
    end
    
    
    if isfield(Info_Correct.MarkersC,'tdtM');
        tdtM = Info_Correct.MarkersC.tdtM(Nbeats,Nchannel)*ParamIn.frequency/1000;
    else
        tdtM = nan(size(X,2),size(X,3));
    end
    
    if isfield(Info_Correct.MarkersC,'tdtm');
        tdtm = Info_Correct.MarkersC.tdtm(Nbeats,Nchannel)*ParamIn.frequency/1000;
    else
        tdtm = nan(size(X,2),size(X,3));
    end
    
    Xiso = nan(size(X,2),size(X,3));%
    Xdtm = nan(size(X,2),size(X,3));
    
    
    
%     if isfield(Info_Correct.MarkersC,'iiTwpos');
%         iiTwpos = Info_Correct.MarkersC.iiTwpos(Nchannel);
%     else
%         iiTwpos= nan(1,size(X,3));
%     end
%     
    tTend_Tneg = nan(size(X,2),size(X,3));
    tTend_Tpos = nan(size(X,2),size(X,3));
    
    if isfield(Info_Correct.MarkersC,'QRSw_fw90');
        QRSw_fw90 = round(Info_Correct.MarkersC.QRSw_fw90(Nbeats,Nchannel)/1000*ParamIn.frequency);
    else
        QRSw_fw90 = nan(size(X,2),size(X,3));
    end
    
    % parameters that maintain the same name
    pp = {'tTpeak','tPw','tTend_defl','tQRSoff','tQRSon','tSw','tQw'};
    for i = 1:length(pp);
        if isfield(Info_Correct.MarkersC,pp{i});
            eval([pp{i},' = round((Info_Correct.MarkersC.(pp{i})(Nbeats,Nchannel) - spikes(Nbeats)*ones(1,length(Nchannel))/ParamIn.frequency*1000)/1000*ParamIn.frequency);']);
        else
            eval([pp{i},' = nan(size(X,2),size(X,3));'])
        end
    end
    
    
    % parameters that do not need any change
    pp = {'Rw_amp','Sw_amp','Qw_amp','QRS_amp','QRS_area','RatioI','Tpeak_amp','ATw','ATw2','Tpeak_amp'};
    for i = 1:length(pp);
        if isfield(Info_Correct.MarkersC,pp{i});
            eval([pp{i},' = Info_Correct.MarkersC.(pp{i})(Nbeats,Nchannel);'])
        else
            eval([pp{i},' = nan(size(X,2),size(X,3));'])
        end
    end
end
% ====== End initialization ================

% ==== manual window
if ~isempty(Info_Correct);
    W_manual = round(Info_Correct.xest/1000*ParamIn.frequency)+DeltaT + [-round(Info_Correct.WinMod/2/1000*ParamIn.frequency):round(Info_Correct.WinMod/2/1000*ParamIn.frequency)];
else
    W_manual = [1:size(X,1)];
end
% ===
if isempty(Info_Correct);
    do_DT = 1;
    do_QRS = 1;
    do_RT_A = 1;
    do_RT_W = 1;
    do_ISO = 1;
    do_Tpeak = 1;
    do_Tend = 1;
    do_P = 1;
else
    do_DT = 0;
    do_QRS = 0;
    do_RT_W = 0;
    do_RT_A = 0;
    do_ISO = 0;
    do_Tpeak = 0;
    do_Tend = 0;
    do_P = 0;
    switch Info_Correct.Marker_name
        case 'AT'
            do_DT = 1;
        case 'RT_Alt'
            do_RT_A = 1;
        case 'RT_Wyatt'
            do_RT_W = 1;
        case  'ISO'
            do_ISO =1;
        case  'Tpeak'
            do_Tpeak =1;
        case 'Tend'
            do_Tend =1;
            do_RT_A =1;
        case 'do_QRS'
            do_QRS = 1;
    end
end

% Hwb = waitbar(0,'Calculating Markers ...');
%% Localize Markers
for i=1:size(X,3) % i=electrode
%     waitbar(i/size(X,3),Hwb);
    if sum(isnan(signals_proc(:,i)))==0
        %         for j = 1:size(X,2) % j=heart beat
        
        
        if do_DT
            Ldt = DTmax+round(80/1000*ParamIn.frequency);
            Wdt = intersect([1:Ldt],W_manual);
            Wdt(Wdt<1 | Wdt>size(X,1)) = [];
            if ~isempty(Wdt)
            ii = Xd2(Wdt,:,i).*Xd2(Wdt+1,:,i)<0 & Xd(Wdt,:,i)<0 & Xd3(Wdt,:,i)>0; % zeros of II derivative & I derivative negative & III derivative > 0
            
            Mask = nan(size(ii));Mask(ii)=1;
            
            [~,kk]= min(Xd(Wdt,:,i).*Mask);
            tm(:,i) = kk+Wdt(1);
            end
        end
        clear ii iim W Ldt Wdt
        
        
        %% QRS
        W = DeltaT+ [-round(QRSw_max/2) : round(QRSw_max/2)];
        W = intersect(W,W_manual);
        W(W<1)=[];
        
        if ~isempty(W)
        iiMax = Xd(W(1:end-1),:,i).*Xd(W(2:end),:,i)<0 & Xd2(W(2:end),:,i)<0; % zeros of I derivative & II derivative negative (V max)
        iiMin = Xd(W(1:end-1),:,i).*Xd(W(2:end),:,i)<0 & Xd2(W(2:end),:,i)>0; % zeros of I derivative & II derivative positive (V min)
        % Polarity of QRS
        if abs(max(nanmean(X(W,:,i),2))) > 0.75*abs(min(nanmean(X(W,:,i),2)))
            QRS_polarity = 1; % R-type
        else
            QRS_polarity = 0; % S-type
        end
        
        % Rwave (max or min within QRS)
        if QRS_polarity
            Maskm = nan(size(iiMax));Maskm(iiMax)=1;
            [m,kk]= max(X(W(2:end),:,i).*Maskm);
            kk(isnan(m)) = nan;
            tdtM(:,i) = kk+W(1);
            x = tdtM(:,i);
            col = [1:size(X,2)];
            XdtM(~isnan(x),i) = diag(X(x(~isnan(x)),col(~isnan(x)),i));
            clear kk Maskm
            
        else
            Maskm = nan(size(iiMin));Maskm(iiMin)=1;
            [m,kk]= min(X(W(2:end),:,i).*Maskm);
            kk(isnan(m)) = nan;
            tdtM(:,i) = kk+W(1);
            x = tdtM(:,i);
            col = [1:size(X,2)];
            XdtM(~isnan(x),i) = diag(X(x(~isnan(x)),col(~isnan(x)),i));
            clear x
        end
        clear m Maskm kk col
        end
        
        
        %% iso time
        if do_ISO
            DTmaxiso = DTmax + round(120/1000*ParamIn.frequency);
            
            Wiso = intersect([1:DTmaxiso],W_manual);
            Wiso(Wiso<1 | Wiso>size(Xd2,1)) = [];
            iidvdtmax = Xd2(Wiso,:,i)<=0 & Xd(Wiso,:,i)>=0; % II derivative negative & I derivative positive
            clear DTmaxiso
            
            Mt = [1:size(iidvdtmax,1)]'*ones(1,size(iidvdtmax,2));
            Mtm = ones(size(Mt,1),1)*tm(:,i)';
            Mtm2 = ones(size(Mt,1),1)*min_RT(:)';
            
            A = Mt>=Mtm & Mt<Mtm2 & Mt>round(60/1000*ParamIn.frequency);
            iidvdtmax=iidvdtmax&A;
            Maskm = nan(size(iidvdtmax));Maskm(iidvdtmax)=1;
            [m,kk] = min(Xd(Wiso,:,i).*Maskm);kk=kk(:);
            kk(isnan(m)) = nan;
            clear m
            
            tdtM(:,i) = kk+W(1);
            x = tdtM(:,i);
            col = [1:size(X,2)];
            XdtM(~isnan(x),i) = diag(X(x(~isnan(x)),col(~isnan(x)),i));
%             XdtM(~isnan(x),i) = diag(X(~isnan(x),[1:size(X,2)],i));
            clear x col
            
            tiso(:,i) = kk + W(1);
            x = tiso(:,i);
            col = [1:size(X,2)];
            Xiso(~isnan(x),i) = diag(X(x(~isnan(x)),col(~isnan(x)),i));
%             Xiso(~isnan(x),i) = diag(X(~isnan(x),[1:size(X,2)],i));
            clear x col
            
            % correction for XdtM: if tdtM==1 & XdtM<Xiso it means that the true XdtM occured before the spike
            XdtM([tdtM(:,i)==1&XdtM(:,i)<Xiso(:,i)],i)=nan;
            
        end
        clear Wiso
        
        %% Repolarization time (Wyatt)
        if do_RT_W||do_RT_A
            clear W
            
            Mt = [1:size(X,1)]'*ones(1,size(X,2));
            MW = Mt <= ones(size(Mt,1),1)*max_RT(:)' &  Mt >= ones(size(Mt,1),1)*min_RT(:)' & Mt>=min(W_manual) & Mt<=max(W_manual);
            
            %             W = intersect(W,W_manual);
            %             W(W>(size(Xd,1)-2))=[];
            
            tt = 1:size(Xd,1)-2;
            ii = Xd2(tt(1:end-1),:,i).*Xd2(tt(2:end),:,i)<0 & Xd(tt(2:end),:,i)>0 & Xd3(tt(2:end),:,i)<=0; % zeros of II derivative & I derivative positive & III derivative < 0(dV/dT max)
            ii_alt = Xd2(tt(1:end-1),:,i).*Xd2(tt(2:end),:,i)<0 & Xd(tt(2:end),:,i)<0 & Xd3(tt(2:end),:,i)>=0;; % zeros of II derivative & I derivative positive & III derivative > 0 (dV/dT min)
            
            ii = ii&MW(tt(1:end-1),:);
            ii_alt = ii_alt&MW(tt(1:end-1),:);
            
            
            if do_RT_W
                if isempty(Info_Correct)
                    M = ones(size(ii));M(Mt(1:size(M,1),:)<ones(size(M,1),1)*tiso(:,i)')=0;
                    ii = ii&M;
                end
                Mask = nan(size(ii));Mask(ii)=1;
                [m,kk] = max(Xd(tt(1:end-1),:,i) .* Mask);
                kk(isnan(m)) = nan;
                tMmax(:,i)=kk+1;
                clear m kk Maks
            end
            
            if do_RT_A
                % downslope
                if isempty(Info_Correct)
                    M = ones(size(ii_alt));M(Mt(1:size(M,1),:)<ones(size(M,1),1)*tMmax(:,i)')=0;
                    ii_alt = ii_alt&M;
                    clear M
                end
                Mask = nan(size(ii));Mask(ii_alt)=1;
                [m,kk] = min(Xd(tt(1:end-1),:,i) .* Mask);
                kk(isnan(m)) = nan;
                tMmin(:,i)=kk+1;
                
                
            end
            
            
            Mt = [1:size(X,1)]'*ones(1,size(X,2));
            
            y = tiso(:,i)';
            y(isnan(y)) = round(170/1000*ParamIn.frequency);
            
            iiA = Mt< ones(size(Mt,1),1)*max([tMmax(:,i),round(tMmin(:,i)*1.2),round(max_RT(:)*.8)],[],2)' & Mt>ones(size(Mt,1),1)*y;
            
            Xbl = X(:,:,i) - ones(size(X,1),1)*Xiso(:,i)';
            Mask = nan(size(iiA));Mask(iiA)=1;
            ATw(:,i) = nansum(Xbl.*Mask)./nansum(abs(Xbl.*Mask));
            
            clear W ii_* ii
        end
        
        if ~isempty(Info_Correct)
            iiTwpos(i) = Info_Correct.MarkersC.iiTwpos(Nchannel(i));
        else
            iiTwpos(i) = nanmedian(ATw(:,i))>=0;
        end
                
        
        %% T-peak
        if do_Tpeak  
            Mt = [1:size(X,1)]'*ones(1,size(X,2));
            MW = Mt <= ones(size(Mt,1),1)*max_RT(:)' &  Mt >= ones(size(Mt,1),1)*min_RT(:)' & Mt>=min(W_manual) & Mt<=max(W_manual);
            
            tt = 1:size(Xd,1)-1;
            iiMax = Xd(tt(1:end-1),:,i).*Xd(tt(2:end),:,i)<0 & Xd2(tt(2:end),:,i)<0; % zeros of I derivative & II derivative negative (V max)
            iiMin = Xd(tt(1:end-1),:,i).*Xd(tt(2:end),:,i)<0 & Xd2(tt(2:end),:,i)>0; % zeros of I derivative & II derivative positive (V min)
            
            iiMax = iiMax&MW(tt(1:end-1),:);
            iiMin = iiMin&MW(tt(1:end-1),:);
            
            
            if iiTwpos(i)
                clear Mask      
                Mask = nan(size(iiMax));Mask(iiMax)=1;
                [m,kk] = max(X(tt(1:end-1),:,i) .* Mask);
                kk(isnan(m)) = nan;
                tTpeak(:,i)=kk+1;
                x = tTpeak(:,i);
                col = [1:size(X,2)];
                Tpeak_amp(~isnan(x),i) = diag(X(x(~isnan(x)),col(~isnan(x)),i));
%                 Tpeak_amp(~isnan(x),i) = diag(X(~isnan(x),[1:size(X,2)],i));
                clear x col Mask
            else
                
                clear M 
                M = ones(size(iiMin));
                M(Mt(1:size(M,1),:)<ones(size(M,1),1)*tiso(:,i)')=0;
                M(Mt(1:size(M,1),:)>ones(size(M,1),1)*tMmax(:,i)')=0;
                iiMin = iiMin&M;
                clear M 
                
                Mask = nan(size(iiMin));Mask(iiMin)=1;
                [m,kk] = min(X(tt(1:end-1),:,i) .* Mask);
                kk(isnan(m)) = nan;
                tTpeak(:,i)=kk+1;
                x = tTpeak(:,i);
                col = [1:size(X,2)];
                Tpeak_amp(~isnan(x),i) = diag(X(x(~isnan(x)),col(~isnan(x)),i));
%                 Tpeak_amp(~isnan(x),i) = diag(X(~isnan(x),[1:size(X,2)],i));
                clear m Mask x
                %
%                 if ~isempty(tt)
%                     [~,i2]=min(X(tt,j,i));
%                     tTpeak(j,i)=tt(i2);
%                     Tpeak_amp(j,i) = X(tTpeak(j,i),j,i);
%                     
%                 end
                
            end
            clear iiMin iiMax W
        end
        
        
        %% Tend deflections [other def]
        
        
    end
    
    
    
    %         end
    
    
    
    %% plot control
    %         figure,plot(X(:,j,i)),hold on,
    %         plot([0 600],[1 1]*prctile(X(1:tiso(j,i),j,i),75),'--k')
    %         plot([0 600],[1 1]*prctile(X(1:tiso(j,i),j,i),25),'--k')
    %         plot([0 600],[1 1]*prctile(X(1:tiso(j,i),j,i),66),'--g')
    %         plot([0 600],[1 1]*prctile(X(1:tiso(j,i),j,i),33),'--g')
    %         plot(tiso(j,i),X(tiso(j,i),j,i),'squarer'),
    %         plot(tMmin(j,i),X(tMmin(j,i),j,i),'or'),
    %         plot(tMmax(j,i),X(tMmax(j,i),j,i),'xr','markersize',10,'linewidth',2),
    %         plot(tTpeak(j,i),X(tTpeak(j,i),j,i),'or','markersize',6,'linewidth',2),
    %         plot(tTe(j,i),X(tTpeak(j,i),j,i),'or','markersize',6,'linewidth',2),
    %         if logical(iiTwpos(i))
    %             plot(tTend_Tpos(j,i),X(round(tTend_Tpos(j,i)),j,i),'squarer','markersize',6,'linewidth',2),
    %         else
    %             plot(tTend_Tneg(j,i),X(round(tTend_Tneg(j,i)),j,i),'squarer','markersize',6,'linewidth',2),
    %         end
    %         plot(tdtM(j,i),X(tdtM(j,i),j,i),'or','markersize',6,'linewidth',2),
    %         plot(tPw(j,i),X(tPw(j,i),j,i),'or','markersize',6,'linewidth',2),
    %
    plotECM = 0;
    
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


RatioI = (Xiso-Xdtm)./(XdtM-Xdtm);


dt = tm+repmat(spikes(Nbeats),[1 size(X,3)]) ; % unwrap
tRw = tdtM+repmat(spikes(Nbeats),[1 size(X,3)])-1 ; % unwrap
tQw = tQw + repmat(spikes(Nbeats),[1 size(X,3)])-1 ; % unwrap
tSw = tSw +repmat(spikes(Nbeats),[1 size(X,3)])-1 ; % unwrap
tQRSon = tQRSon + repmat(spikes(Nbeats),[1 size(X,3)])-1 ; % unwrap
tQRSoff = tQRSoff +repmat(spikes(Nbeats),[1 size(X,3)])-1 ; % unwrap
rt_up = tMmax+repmat(spikes(Nbeats),[1 size(X,3)]) ; % Unwrap
rt_down = tMmin+repmat(spikes(Nbeats),[1 size(X,3)]) ; % Unwrap
isot = tiso+repmat(spikes(Nbeats),[1 size(X,3)]) ; % unwrap
tTpeak = tTpeak+repmat(spikes(Nbeats),[1 size(X,3)])-1; % unwrap
tTend_Tpos = tTend_Tpos+repmat(spikes(Nbeats),[1 size(X,3)]) ; % unwrap
tTend_Tneg = tTend_Tneg+repmat(spikes(Nbeats),[1 size(X,3)]) ; % unwrap
tTend_defl = tTend_defl+repmat(spikes(Nbeats),[1 size(X,3)])-1 ; % unwrap
tPw = tPw+repmat(spikes(Nbeats),[1 size(X,3)]) - 1;
tPwave_int = tdtM - tPw;


%% ---
do_plot = 0;
if do_plot
    i = 1;
    figure,
    plot(signals_proc(:,i)),hold on,
    plot(tRw(~isnan(tRw(:,i)),i),signals_proc(tRw(~isnan(tRw(:,i)),i),i),'or')
    ,plot(tSw(~isnan(tSw(:,i)),i),signals_proc(tSw(~isnan(tSw(:,i)),i),i),'squareg'),
    plot(tQw(~isnan(tQw(:,i)),i),signals_proc(tQw(~isnan(tQw(:,i)),i),i),'squareg')
    plot(tTpeak(~isnan(tTpeak(:,i)),i),signals_proc(tTpeak(~isnan(tTpeak(:,i)),i),i),'*k')
    plot(tPw(~isnan(tPw(:,i)),i),signals_proc(tPw(~isnan(tPw(:,i)),i),i),'om')
    if iiTwpos(i)
        xx = round(tTend_Tpos(~isnan(tTend_Tpos(:,i)),i));xx(xx>size(signals_proc,1))=[];
        plot(xx,signals_proc(xx,i),'squarek')
    else
        xx = round(tTend_Tneg(~isnan(tTend_Tneg(:,i)),i));xx(xx>size(signals_proc,1))=[];
        plot(xx,signals_proc(xx,i),'squarek')
    end
    plot(tTend_defl(~isnan(tTend_defl(:,i)),i),signals_proc(round(tTend_defl(~isnan(tTend_defl(:,i)),i)),i),'^k')
    
    set(gca,'xtick',spikes+DeltaT,'xgrid','on')
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
    'tRw [ms]: R-wave (defined as the prominent wave independently of polarity)',...,
    'tSw [ms]: S-wave (defined as the wave following the prominent one)',...,
    'tQw [ms]: Q-wave (defined as the wave preceding the prominent one)',...,
    'tPwave_int [ms]: PR interval',...,
    'tPw [ms]: Instant of Pwave peak',...,
    'RatioI: [V(isot)-V(tdtm)]/[V(tdtM)-V(tdtm)]',...,
    'ParamIn: parameters of analysis'};


iiTwpos = logical(iiTwpos); % true if T-wave is positive, false otherwise
% from samples to milli-seconds
Markers.tTpeak = tTpeak/ParamIn.frequency*1000;

Markers.tPwave_int = tPwave_int/ParamIn.frequency*1000;
Markers.tPw = tPw/ParamIn.frequency*1000;

Markers.tTend = nan(size(dt));
Markers.tTend(:,iiTwpos) = tTend_Tpos(:,iiTwpos)/ParamIn.frequency*1000;
Markers.tTend(:,~iiTwpos) = tTend_Tneg(:,~iiTwpos)/ParamIn.frequency*1000;
Markers.tTend_defl = tTend_defl/ParamIn.frequency*1000;
% Markers.tTend(end,[Markers.tTend(end,:)>size(signals_proc,1)]) = nan;
Markers.tTend(Markers.tTend/1000*ParamIn.frequency>size(signals_proc,1)) = nan;
Markers.tTend(Markers.tTend<0)=nan;

Markers.rt_Wyatt = rt_up/ParamIn.frequency*1000;
Markers.rt_Alternative(:,~iiTwpos) =rt_up(:,~iiTwpos)/ParamIn.frequency*1000;
Markers.rt_Alternative(:,iiTwpos) = rt_down(:,iiTwpos)/ParamIn.frequency*1000;

Markers.tSw = tSw/ParamIn.frequency*1000;
Markers.tRw = tRw/ParamIn.frequency*1000;
Markers.tQw = tQw/ParamIn.frequency*1000;
Markers.tQRSon = tQRSon/ParamIn.frequency*1000;
Markers.tQRSoff = tQRSoff/ParamIn.frequency*1000;
Markers.QRSw_fw90 = QRSw_fw90/ParamIn.frequency*1000;
Markers.QRSw = Markers.tSw-Markers.tQw;
Markers.tQRSon = tQRSon/ParamIn.frequency*1000;
Markers.tQRSoff = tQRSoff/ParamIn.frequency*1000;


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
Markers.m_function = mfilename('fullpath');
Markers.Legend = leg;
Markers.Tpeak_amp = Tpeak_amp;
Markers.Rw_amp = Rw_amp;
Markers.Sw_amp = Sw_amp;
Markers.Qw_amp = Qw_amp;
Markers.QRS_amp = QRS_amp; % peak to peak (range(QRS))
Markers.QRS_area = QRS_area; % Area from Qw to Sw



% Cancel possible negative values in first beat
vv = fieldnames(Markers);
for i = 1:length(vv)
    if size(Markers.(vv{i}),1)==length(spikes)
        if isnumeric(Markers.(vv{i}))
            Markers.(vv{i})(1,(Markers.(vv{i})(1,:))<1)=nan;
        end
    end
end


% close(Hwb)



% Ex DTRT_mo_gui_bio_bank_v3.m
function [Markers,Hwb] = DTRT_mo_gui_v2(signals_proc,spikes,ParamIn,Info_Correct);
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
L=min([max(diff(spikes)),maxCL]) + Hend;
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
    
    
    
    if isfield(Info_Correct.MarkersC,'iiTwpos');
        iiTwpos = Info_Correct.MarkersC.iiTwpos(Nchannel);
    else
        iiTwpos= nan(1,size(X,3));
    end
    
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

Hwb = waitbar(0,'Calculating Markers ...');
%% Localize Markers
for i=1:size(X,3) % i=electrode
    waitbar(i/size(X,3),Hwb);
    if sum(isnan(signals_proc(:,i)))==0
        Ldt = DTmax+round(80/1000*ParamIn.frequency);
        
        Wdt = intersect([1:Ldt],W_manual);
        Wdt(Wdt<1 | Wdt>size(X,1)) = [];
        
        ii = Xd2(Wdt,:,i).*Xd2(Wdt+1,:,i)<0 & Xd(Wdt,:,i)<0 & Xd3(Wdt,:,i)>0; % zeros of II derivative & I derivative negative & III derivative > 0
        iim = Xd(Wdt,:,i).*Xd(Wdt+1,:,i)<0; % zeros of I derivative
        if do_DT
            %% Depolarization time
            if sum(ii(:))>0
                for j = 1:size(ii,2) % j=heart beat
                    tt = Wdt(ii(:,j));
                    if ~isempty(tt)
                        [~,kk]=min(Xd(tt,j,i));
                        tm(j,i) = tt(kk);
                        ttM = Wdt(iim(:,j));
                        it1 = find(ttM<tt(kk) & Xd2(ttM,j,i)'<0);
                        it2 = find(ttM>tt(kk) & ttM<tt(kk)+100 & Xd2(ttM,j,i)'>0);
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
                    
                end
            end
        end
        clear ii iim W Ldt Wdt
        
        if do_QRS
            %%
            W = DeltaT+ [-round(QRSw_max/2) : round(QRSw_max/2)];
            W = intersect(W,W_manual);
            W(W<1)=[];
            
            iiMax = Xd(W(1:end-1),:,i).*Xd(W(2:end),:,i)<0 & Xd2(W(2:end),:,i)<0; % zeros of I derivative & II derivative negative (V max)
            iiMin = Xd(W(1:end-1),:,i).*Xd(W(2:end),:,i)<0 & Xd2(W(2:end),:,i)>0; % zeros of I derivative & II derivative positive (V min)
            % Polarity of QRS
            if abs(max(nanmean(X(W,:,i),2))) > 0.75*abs(min(nanmean(X(W,:,i),2)))
                QRS_polarity = 1; % R-type
            else
                QRS_polarity = 0; % S-type
            end
            
            for j = 1:size(X,2) % j=heart beat
                
                %% QRS analysis
                % Rwave
                if QRS_polarity
                    tt = W(iiMax(:,j));
                    if ~isempty(tt)
                        [~,kk]=max(X(tt,j,i));
                        tdtM(j,i) = tt(kk) + 1;
                        XdtM(j,i) = X(tdtM(j,i),j,i);
                    end
                else
                    tt = W(iiMin(:,j));
                    if ~isempty(tt)
                        [~,kk]=min(X(tt,j,i));
                        
                        tdtM(j,i) = tt(kk) + 1;
                        XdtM(j,i) = X(tdtM(j,i),j,i);
                    end
                end
                
                clear tt
                % QRS width (robust)
                tp1 = [];tp2=[];
                if ~isnan(tdtM(j,i))
                    HH = tdtM(j,i)+[round(-QRSw_max/2+1) : round(QRSw_max/2)];
                    HH(HH<1)=[];
                    if QRS_polarity
                        fw90m = X(HH,j,i) > 0.10*XdtM(j,i);
                        tp1 = [1 HH(find(diff(fw90m)==1))];tp1(tp1>tdtM(j,i))=[];
                        tp2 = [HH(find(diff(fw90m)==-1)) HH(end)];tp2(tp2<tdtM(j,i))=[];
                    else
                        fw90m = X(HH,j,i) < 0.10*XdtM(j,i);
                        tp1 = [1 HH(find(diff(fw90m)==1))];tp1(tp1>tdtM(j,i))=[];
                        tp2 = [HH(find(diff(fw90m)==-1)) HH(end)];tp2(tp2<tdtM(j,i))=[];
                    end
                    [~,i1] = min(abs(tp1-tdtM(j,i)));tp1 = tp1(i1);
                    [~,i2] = min(abs(tp2-tdtM(j,i)));tp2 = tp2(i2);
                    clear i1 i2
                    QRSw_fw90(j,i) = tp2 - tp1;  % samples
                    clear HH
                end
                tp12 = round((tp1+tp2)/2);
                % Q-wave & T-on
                if ~isempty(tp12)
                    HH =  tp12 + [-round(90/1000*ParamIn.frequency) : round(10/1000*ParamIn.frequency)];
                    HH(HH<1)=[];
                    if QRS_polarity
                        tt = find(Xd(HH(2:end),j,i).*Xd(HH(1:end-1),j,i)<0 & Xd2(HH(1:end-1),j,i)>=0);
                        tton = find(Xd(HH(2:end),j,i).*Xd(HH(1:end-1),j,i)<=0 & Xd2(HH(1:end-1),j,i)<=0);
                        
                        xx = find(Xd2(HH(1:end-1),j,i).*Xd2(HH(2:end),j,i) <=0 & Xd(HH(1:end-1),j,i)<0);
                        [m,im] = max( Xd(HH(xx),j,i));
                        tton2 = HH(xx(im));
                    else
                        tt = find(Xd(HH(2:end),j,i).*Xd(HH(1:end-1),j,i)<0 & Xd2(HH(1:end-1),j,i)<=0);
                        tton = find(Xd(HH(2:end),j,i).*Xd(HH(1:end-1),j,i)<0 & Xd2(HH(1:end-1),j,i)>=0);
                        xx = find(Xd2(HH(1:end-1),j,i).*Xd2(HH(2:end),j,i) <=0 & Xd(HH(1:end-1),j,i)>0);
                        [m,im] = min( Xd(HH(xx),j,i));
                        tton2 = HH(xx(im));
                        
                    end
                    if ~isempty(tt)
                        % Qw
                        [~,iq] = min(abs(HH(tt)-tp12));tt = tt(iq);
                        % Ton
                        tton(tton>tt)=[];
                        [~,iq] = min(abs(tton-tt));tton = tton(iq);
                        if ~isempty(tt)
                            tQw(j,i) = tt+HH(1);
                        end
                        if ~isempty(tton)
                            tQRSon(j,i) = tton+HH(1);
                        else
                            if m/max(abs(Xd(HH,j,i)))<0.01
                                tQRSon(j,i) = tton2;
                            end
                        end
                    end
                end
                
                clear tt HH tton iq
                
                % S-wave % QRS_off
                if ~isempty(tp12)
                    HH =  tp12 + [-round(-10/1000*ParamIn.frequency) : round(90/1000*ParamIn.frequency)];
                    HH(HH<1)=[];
                    if QRS_polarity
                        tt = find(Xd(HH(2:end),j,i).*Xd(HH(1:end-1),j,i)<0 & Xd2(HH(1:end-1),j,i)>=0);
                        ttoff = find(Xd(HH(2:end),j,i).*Xd(HH(1:end-1),j,i)<0 & Xd2(HH(1:end-1),j,i)<=0);
                        
                        xx = find(Xd2(HH(1:end-1),j,i).*Xd2(HH(2:end),j,i) <=0 & Xd(HH(1:end-1),j,i)>0);
                        [m,im] = min( Xd(HH(xx),j,i));
                        ttoff2 = HH(xx(im));
                        
                    else
                        tt = find(Xd(HH(2:end),j,i).*Xd(HH(1:end-1),j,i)<0 & Xd2(HH(1:end-1),j,i)<=0);
                        ttoff = find(Xd(HH(2:end),j,i).*Xd(HH(1:end-1),j,i)<0 & Xd2(HH(1:end-1),j,i)>=0);
                        
                        xx = find(Xd2(HH(1:end-1),j,i).*Xd2(HH(2:end),j,i) <=0 & Xd(HH(1:end-1),j,i)<0);
                        [m,im] = max( Xd(HH(xx),j,i));
                        ttoff2 = HH(xx(im));
                        
                    end
                    if ~isempty(tt)
                        % Sw
                        [~,is] = min(abs(HH(tt)-tp12));tt = tt(is);
                        % tQRSoff
                        ttoff(ttoff<tt)=[];
                        [~,is] = min(abs(ttoff-tt));ttoff = ttoff(is);
                        
                        if ~isempty(tt)
                            tSw(j,i) = tt+HH(1);
                        end
                        if ~isempty(ttoff)
                            tQRSoff(j,i) = ttoff+HH(1);
                        else
                            
                            if m/max(abs(Xd(HH,j,i)))<0.01
                                tQRSoff(j,i) = ttoff2;
                            end
                        end
                    end
                end
                
                clear tt HH is
                %             % plot to check
                %             figure(1),
                %             hold off
                %             plot(X(:,j,i)),hold on,
                %             plot(W,X(W,j,i),'g'),plot(tdtM(j,i),XdtM(j,i),'squarer'),
                %             plot(tQw(j,i),X(tQw(j,i),j,i),'or'),plot(tSw(j,i),X(tSw(j,i),j,i),'or')
                %             plot(tQRSon(j,i),X(tQRSon(j,i),j,i),'xr'),plot(tQRSoff(j,i),X(tQRSoff(j,i),j,i),'xr')
                %             title(['chan=',num2str(j)])
                %             %             pause
                
                % - QRS-params
                if ~isnan(XdtM(j,i))
                    Rw_amp(j,i) = XdtM(j,i);
                end
                if ~isnan(tSw(j,i))
                    Sw_amp(j,i) = X(tSw(j,i),j,i);
                end
                if ~isnan(tQw(j,i))
                    Qw_amp(j,i) = X(tQw(j,i),j,i);
                end
                t1 = min([tQRSon(j,i),tQw(j,i)]);t2 = max([tQRSoff(j,i),tSw(j,i)]);
                if ~isnan(t1+t2)
                    QRS_amp(j,i) = range(X(t1:t2,j,i));
                end
                clear t1 t2
                if ~isnan(tQw(j,i))&~isnan(tSw(j,i))
                    QRS_area(j,i) = nanmean(abs( X(tQw(j,i):tSw(j,i),j,i) - X(tQw(j,i),j,i)));
                end
            end
        end
        clear ii iim Ldt W
        
        %% iso time
        if do_ISO
            DTmaxiso = DTmax + round(120/1000*ParamIn.frequency);
            
            Wiso = intersect([1:DTmaxiso],W_manual);
            Wiso(Wiso<1 | Wiso>size(Xd2,1)) = [];
            iidvdtmax = Xd2(Wiso,:,i)<=0 & Xd(Wiso,:,i)>=0; % II derivative negative & I derivative positive
            clear DTmaxiso
            if sum(iidvdtmax(:))>0
                for j = 1:size(iidvdtmax,2) % j=heart beat
                    ttdvmax = Wiso(iidvdtmax(:,j))+1;
                    if ~isempty(ttdvmax)
                        ttdvmax(ttdvmax<tm(j,i)|ttdvmax<tdtm(j,i)+15|ttdvmax < round(80/1000*ParamIn.frequency)| ttdvmax > min_RT(j))=[];
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
        clear Wiso
        
        %%
        %% Repolarization time (Wyatt)
        if do_RT_W||do_RT_A
            for j=1:size(X,2) % beats
                
                W = min_RT(j) : max_RT(j); % 03/11/2017 inside the loop because in theory this can change for every beat
                W = intersect(W,W_manual);
                W(W>(size(Xd,1)-2))=[];
                
                ii = Xd2(W(1:end-1),j,i).*Xd2(W(2:end),j,i)<0 & Xd(W(2:end),j,i)>0 & Xd3(W(2:end),j,i)<=0; % zeros of II derivative & I derivative positive & III derivative < 0(dV/dT max)
                ii_alt = Xd2(W(1:end-1),j,i).*Xd2(W(2:end),j,i)<0 & Xd(W(2:end),j,i)<0 & Xd3(W(2:end),j,i)>=0;; % zeros of II derivative & I derivative positive & III derivative > 0 (dV/dT min)
                iiMin = Xd(W(1:end-1),j,i).*Xd(W(2:end),j,i)<0 & Xd2(W(2:end),j,i)>0; % zeros of I derivative & II derivative positive (V min)
                iiMax = Xd(W(1:end-1),j,i).*Xd(W(2:end),j,i)<0 & Xd2(W(2:end),j,i)<0; % zeros of I derivative & II derivative negative (V max)
                if do_RT_W
                    
                    tt = W(ii(:));
                    % not too close to the spikes
                    if squeeze(nanmedian(sum(X(W,:,i)),2))>0 % Only for positive T-wave
                        if j<size(X,2)
                            tt(tt>spikes(j+1)-spikes(j)- round(10/1000*ParamIn.frequency))=[]; % Upstroke no too close to spike
                        end
                    end
                    if isempty(Info_Correct)
                        tt(tt<=tiso(j,i))=[];
                    end
                    %                     tt(tt>tm(j,i)+max_ARI|tt<tm(j,i)+min_ARI) = []; % control over ARI
                    
                    if ~isempty(tt)
                        if squeeze(nanmedian(sum(X(W,:,i)),2))>0 % Only for positive T-wave
                            if ~isempty(max(find([Xd2(tt,j,i)-Xd2(tt-1,j,i)]>0)))&max(find([Xd2(tt,j,i)-Xd2(tt-1,j,i)]>0))<length(tt)& max(find([Xd2(tt,j,i)-Xd2(tt-1,j,i)]>0))>1
                                tt(1:max(find([Xd2(tt,j,i)-Xd2(tt-1,j,i)]>0)))=[]; % eliminate fiducial points which occurr at same time or before the first change from convex to concave (can be part of a slow depolarization)
                            end
                        end
                        
                        [~,ii2] = max(Xd(tt,j,i));
                        tMmax(j,i)=tt(ii2);
                    end
                end
                if do_RT_A
                    % downslope
                    ttalt = W(ii_alt(:));
                    if isempty(Info_Correct)
                        if j<size(X,2)
                            ttalt(ttalt<tMmax(j,i) | spikes(j+1)-spikes(j)-ttalt < round(10/1000*ParamIn.frequency))=[];
                        else
                            ttalt(ttalt<tMmax(j,i))=[];
                        end
                    end
                    if ~isempty(ttalt)
                        [p1down,ii2] = min(Xd(ttalt,j,i));
                        tMmin(j,i)=ttalt(ii2);
                    end
                end
                
                
                % Area Twave
                tr = max([tMmax(j,i),round(tMmin(j,i)*1.2),round(max_RT(j)*.8)]);
                tr = min([tr,size(X,1)]);
                if ~isnan(tiso(j,i)+tr); % mo 17/06/2015
                    Xbl = X(:,j,i)-X(tiso(j,i),j,i);
                    ATw(j,i) = sum(Xbl(tiso(j,i):tr))/nansum(abs(Xbl(tiso(j,i):tr)));
                else
                    xx = round(prctile(tm(j,:),95)+round(50/1000*ParamIn.frequency)) : round(prctile(tMmin(j,:),95)+round(50/1000*ParamIn.frequency));
                    if isnan(xx)
                        xx = round(round(170/1000*ParamIn.frequency)) : round(400/1000*ParamIn.frequency);
                        xx(xx>size(X,1)) = [];
                    end
                    A = detrend(X(xx,j,i));
                    A = A-nanmean( A(1:round(50/1000*ParamIn.frequency) ));
                    ATw(j,i) = nansum(A)/nansum(abs(A));
                end
                
                
                % Area Twave2
                if ~isnan(tiso(j,i));
                    A = X(tiso(j,i):end,j,i);A(isnan(A))=[];
                    Xbl = detrend(A);clear A
                    Xbl=Xbl-Xbl(1);
                    ATw2(j,i) = sum(Xbl)/nansum(abs(Xbl));
                end
                
            end
            clear W ii_* ii
        end
        iiTwpos(i) = nanmedian(ATw(:,i))>=0;
        
        %% T-end
        if do_Tend
            for j=1:size(X,2)
                W = min_RT(j) : max_RT(j); %
                W = intersect(W,W_manual);
                W(W>(size(Xd,1)-2))=[];
                
                iiMin = Xd(W(1:end-1),j,i).*Xd(W(2:end),j,i)<0 & Xd2(W(2:end),j,i)>0; % zeros of I derivative & II derivative positive (V min)
                iiMax = Xd(W(1:end-1),j,i).*Xd(W(2:end),j,i)<0 & Xd2(W(2:end),j,i)<0; % zeros of I derivative & II derivative negative (V max)
                
                
                % positive Twave
                if iiTwpos(i) && ~isnan(tMmin(j,i))
                    x0 = tMmin(j,i); % samples
                    p1down = Xd(x0,j,i);
                    y0 = X(tMmin(j,i),j,i);
                    %                     tte = W(find(iiMin(:,j)))+1;
                    tte = min(W((iiMin(:)))+1); % 25/11/2016
                    tte(tte<x0)=[];
                    if j==length(spikes)
                        tte2 = round(min(spikes(j)+nanmedian(spikes)-200/1000*ParamIn.frequency -spikes(j),size(Xd2,1)));
                    else
                        %                         tte2 = spikes(j+1)- round(200/1000*ParamIn.frequency) -spikes(j);
                        tte2 = min(spikes(j+1)-round(200/1000*ParamIn.frequency)-spikes(j),size(Xd3,1));
                    end
                    %             if ~isnan(tte)
                    %             tTend_defl(j,i) = tte(1);
                    %             end
                    Tm = tte : tte2;
                    Tm(abs(Xd(tte:tte2,j,i))>abs(0.1*p1down))=[];
                    it = find(abs(diff(Tm))>1);
                    if ~isempty(it)
                        Tm(it(1)+1:end)=[];
                    end
                    V0 = nanmean(X(Tm,j,i));
                    if isnan(V0); V0=0; end
                    clear it Tm tte tte2
                    if ~isnan(x0)
                        %                     rV0 = roots([p1down  (y0-p1down*x0-V0)]);
                        rV0 = -(y0-p1down*x0-V0)/p1down;
                        if rV0<size(X,1) % no longer than cycle length
                            %                             tTend_Tpos(j,i) = rV0/ParamIn.frequency*1000; % time for which the line tangent to tMin is equal to zero (for positive T waves)
                            tTend_Tpos(j,i) = round(rV0); % time for which the line tangent to tMin is equal to zero (for positive T waves)
                        end
                    end
                end
                %             % plot
                %                 t = [0:size(X,1)-1];
                %                 y = p1down*(t-x0)+y0;
                %                 figure,
                %                 plot(t,X(:,j,i)),hold on,
                %                 plot(t,y,'r');
                %                 plot(round(rV0)*[1 1],get(gca,'ylim'),'--r')
                %                 grid on
                %                 clear x0 y0 tte* Tm
                
                % negative Twaves
                if ~iiTwpos(i) && ~isnan(tMmax(j,i))
                    
                    x0 = tMmax(j,i);
                    p1up = Xd(x0,j,i);
                    y0 = X(tMmax(j,i),j,i);
                    tte = min(W((iiMax(:)))+1);
                    tte(tte<x0)=[];
                    if j==length(spikes)
                        tte2 = round(min(spikes(j)+nanmedian(spikes)-200/1000*ParamIn.frequency -spikes(j),size(Xd2,1)));
                    else
                        tte2 = min(spikes(j+1) - round(200/1000*ParamIn.frequency)-spikes(j),size(Xd3,1));
                    end
                    %             if ~isnan(tte)
                    %             tTend_defl(j,i) = tte(1);
                    %             end
                    Tm = tte:tte2;
                    Tm(abs(Xd(tte:tte2,j,i))>abs(0.1*p1up))=[];
                    it = find(abs(diff(Tm))>1);
                    if ~isempty(it)
                        Tm(it(1)+1:end)=[];
                    end
                    V0 = nanmean(X(Tm,j,i));
                    if isnan(V0); V0=0; end
                    clear it Tm tte tte2
                    if ~isnan(x0)
                        %                     rV0 = roots([p1up (y0-p1up*x0-V0)]);
                        rV0 = -(y0-p1up*x0-V0)/p1up;
                        if rV0<size(X,1) % no longer than cycle length
                            %                         tTend_Tneg(j,i) = rV0/ParamIn.frequency*1000; % time for which the line tangent to tMin is equal to zero (for positive T waves)
                            tTend_Tneg(j,i) = round(rV0); % time for which the line tangent to tMin is equal to zero (for positive T waves)
                            
                        end
                    end
                    %                                      % plot
                    %                             t = [0:size(X,1)-1];
                    %                             y = p1up*(t-x0)+y0;
                    %                             figure,
                    %                             plot(t,X(:,j,i)),hold on,
                    %                             plot(t,y,'r');
                    %                 plot(round(rV0)*[1 1],get(gca,'ylim'),'--r')
                    %                             grid on
                    
                    
                    %             tt = find(ii(:,j))+W(1)-1;
                    %             if ~isempty(tt)
                    %                 [p1up,ii2] = max(Xd(tt,j,i));
                    %             end
                    %
                    %             if ~isnan(tMmax(j,i))/ParamIn.frequency*1000
                    %                 r = roots([p1up  (X(round(tMmax(j,i)/1000*ParamIn.frequency),j,i)-p1up*tMmax(j,i))]);
                    %                 if r<size(X,1) % no longer than cycle length
                    %                     tTend_Tneg(j,i) = r; % time for which the line tangent to tMax is equal to zero (for positive T waves)
                    %                 end
                    %             end
                end
            end
            %% %%%%%%%%%%%%%%%%%%
            
            
            %%
        end
        
        
        
        %% T-peak
        if do_Tpeak
            for j=1:size(X,2)
                
                W = min_RT(j) : max_RT(j); % 10/05/2017
                W = intersect(W,W_manual);
                W(W>(size(Xd,1)-1))=[];
                
                iiMax = Xd(W(1:end-1),j,i).*Xd(W(2:end),j,i)<0 & Xd2(W(2:end),j,i)<0; % zeros of I derivative & II derivative negative (V max)
                iiMin = Xd(W(1:end-1),j,i).*Xd(W(2:end),j,i)<0 & Xd2(W(2:end),j,i)>0; % zeros of I derivative & II derivative positive (V min)
                
                if iiTwpos(i)
                    tt = W((iiMax(:)))+1;
%                     tt(tt<=tMmax(j,i)|tt>=tMmin(j,i))=[]; Within RT_Wyatt
%                     and RT_Alternative
                    if ~isempty(tt)
                        [~,i2]=max(X(tt,j,i));
                        tTpeak(j,i) = tt(i2);
                        Tpeak_amp(j,i) = X(tTpeak(j,i),j,i);
                    end
                    clear tt i2
                    %                     tte = Wte((iiMin_te(:,j)))+1;
                    %                     tte(tte<tTpeak(j,i))=[];
                    %                     if ~isnan(tte)
                    %                         tTend_defl(j,i) = tte(1);
                    %                     end
                    %                     clear tte
                else% negative Twave
                    tt = W((iiMin(:)))+1;
                    tt(tt>=tMmax(j,i)|tt<tiso(j,i))=[];
                    if ~isempty(tt)
                        [~,i2]=min(X(tt,j,i));
                        tTpeak(j,i)=tt(i2);
                        Tpeak_amp(j,i) = X(tTpeak(j,i),j,i);
                        %                         clear tt i2
                        %                         tte = Wte((iiMax_te(:,j)))+1;
                        %                         tte(tte<tTpeak(j,i))=[];
                        %                         if ~isnan(tte)
                        %                             tTend_defl(j,i) = tte(1);
                        %                         end
                        %                         clear tte
                    end
                    
                end
            end
            clear iiMin iiMax W
        end
        
        
        %% Tend deflections [other def]
        if do_Tend
            for j=1:size(X,2)
                %             W = [min(min_RT) : max(max_RT)+ round(0.050*ParamIn.frequency)]; % windows in which rep time can be estimated
                W = [min_RT(j) : max_RT(j)+round(0.050*ParamIn.frequency)]; % 10/05/2017
                W(W>size(Xd3,1)) = [];
                iiD1 = (Xd(W(1:end-1),:,i)).*(Xd(W(2:end),:,i))<0 ; % zeros of I derivative
                %             te = W(iiD1(:,j))+1;
                %             y = detrend(X(W,j,i));
                %             [~,iim] = max(abs(y(te-W(1))));
                %             atp = X(te(iim),j,i);
                %             tp=te(iim);
                tp = tTpeak(j,i);
                if ~isnan(tp)
                    atp = X(tp,j,i);
                else
                    atp = max(abs(X(:,j,i)));
                end
                iiD2 = [Xd2(W(1:end-1),:,i).*Xd2(W(2:end),:,i)<0 & Xd(W(2:end),:,i)>=0];
                
                if atp<=0 % positive Twave = discard positive peaks
                    iiD1 =[iiD1 & Xd2(W(2:end),:,i)<=0];
                else
                    iiD1 =[iiD1 & Xd2(W(2:end),:,i)>=0];
                end
                
                tmD1 = W(iiD1(:,j))+1;
                if ~isempty(tp)
                    tmD1(tmD1<=tp+0.020*ParamIn.frequency)=[]; % tPeak
                    tmD1( abs(X(tmD1,j,i)) > 0.3*abs(atp))=[]; % Amplitude
                end
                if ~isempty(tmD1)
                    tmD1 = tmD1(1);
                end
                
                
                tmD2 = W(iiD2(:,j))+1;
                if ~isempty(tp)
                    tmD2(tmD2<=tp+0.020*ParamIn.frequency)=[]; % tPeak
                    tmD2( abs(X(tmD2,j,i)) > 0.3*abs(atp))=[]; % Amplitude
                    
                end
                if ~isempty(tmD1)
                    tmD2(tmD2>=tmD1)=[]; %
                end
                
                te = sort([tmD1 tmD2],'ascend');
                
                tmm = nan(size(te));
                for kk = 1:length(te)
                    tmm(kk) = nanmean(abs(Xd(te(kk)-round(0.005*ParamIn.frequency) : min([size(Xd,1),te(kk)+round(0.005*ParamIn.frequency)]),j,i)));
                end
                [m,iim] = min(tmm);
                
                te =te(iim);
                if ~isempty(te)
                    tTend_defl(j,i) = te(1);
                end
            end
        end
        clear W tmm iim tmD2 tmD1 atp y te
        %     figure,
        %     plot(X(:,1,i)),hold on,
        %     plot(tmD1,X(tmD1,1,i),'ok','markersize',8)
        %     plot(tmD2,X(tmD2,1,i),'or','markersize',8)
        %     plot(tm,X(tm,1,i),'xr','markersize',14)
        %     xlim([500 1600])
        %     pause
        
        
        %% P-wave
        clear W
        if do_P && DeltaT>=120/1000*ParamIn.frequency % samples
            
            W = [1 : DeltaT]; % windows in which rep time can be estimated
            W = intersect(W,W_manual);
            
            iiMax = Xd(W(1:end-1),:,i).*Xd(W(2:end),:,i)<0 & Xd2(W(2:end),:,i)<0; % zeros of I derivative & II derivative negative (V max)
            iiMin = Xd(W(1:end-1),:,i).*Xd(W(2:end),:,i)<0 & Xd2(W(2:end),:,i)>0; % zeros of I derivative & II derivative positive (V min)
            for j=1:size(ii,2)
                if iiTwpos(i) % we assume that P and T wave has the same polarity (it seems to be right in this configuration)
                    tt = W(find(iiMax(:,j)))+1;
                    tt(tt>=tdtM(j,i))=[];
                    tt( (tdtM(j,i)-tt)<tPi_min | (tdtM(j,i)-tt)>tPi_max)=[];
                    if ~isempty(tt)
                        [~,i2]=max(X(tt,j,i));
                        tPw(j,i)=tt(i2);
                    end
                    clear tt i2
                else% negative Twave
                    tt = W(find(iiMin(:,j)))+1;
                    tt(tt>=tdtM(j,i))=[];
                    tt( (tdtM(j,i)-tt)<tPi_min | (tdtM(j,i)-tt)>tPi_max)=[];
                    if ~isempty(tt)
                        [~,i2]=min(X(tt,j,i));
                        tPw(j,i)=tt(i2);
                        clear tt i2
                    end
                    
                end
            end
            clear W j
            
        end
        
        
        
    end
    
    
    
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


close(Hwb)



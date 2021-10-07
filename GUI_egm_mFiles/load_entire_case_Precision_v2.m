clear
close all
addpath 'E:\UCL\Scripts_all\Scripts_mo\GUI_egm\GUI_egm_mFiles'
%% Load mesh
[FILENAME, PATHNAME] = uigetfile('*DxLandmarkGeo.xml','Select MESH'); % select file
meshname = [PATHNAME,FILENAME];
[ParamOut] = load_DxLandmarkGeo_Precision(meshname,1);

% -
X = dir([PATHNAME,'DxL*.csv']);
S = [];
labels = {};
xyz = []; xyz_surfP = [];
S_ref = [];S_spare1 = [];
S_spare2 = [];S_spare3 = [];
h = waitbar(0,'Converting data ...');
for ip = 1:length(X);
    waitbar(ip/length(X),h);
    filename = [PATHNAME,X(ip).name];
    [s,p] = load_DxL_Precision(filename);
    
    S = [S s];
    labels = cat(2,labels,p.Label);
    xyz = [xyz p.xyz'];
    xyz_surfP = [xyz_surfP p.xyz_surfP'];
    
    S_ref = [S_ref p.Other_signals.signals_ref];
    S_spare1 = [S_spare1 p.Other_signals.signals_spare1];
    S_spare2 = [S_spare2 p.Other_signals.signals_spare2];
    S_spare3 = [S_spare3 p.Other_signals.signals_spare3];
    
end
close(h);
ParamAnalysis.N_ori = size(S,2);
ParamSig.frequency = p.Frequency;

% = Eliminating "zeros"
ParamSig.Label = labels; clear labels
iiko = mean(S==0)>0.99;
ParamSig.Label(iiko) = [];
S(:,iiko) = [];
xyz(:,iiko) = [];
xyz_surfP(:,iiko) = [];
S_ref(:,iiko) = [];
S_spare1(:,iiko) = [];
S_spare2(:,iiko) = [];
S_spare3(:,iiko) = [];
ParamAnalysis.Removed_equal_to_zero = sum(iiko);
clear iiko

% removing duplicates
[A,iA] = unique(ParamSig.Label,'stable');
ParamAnalysis.Removed_duplicates = size(S,2)-length(iA);
ParamSig.Label = A; clear A
signals = S(:,iA);
xyz = xyz(:,iA);
xyz_surfP =  xyz_surfP(:,iA);
S_ref = S_ref(:,iA);
S_spare1 = S_spare1(:,iA);
S_spare2 = S_spare2(:,iA);
S_spare3 = S_spare3(:,iA);

othersignals.S_spare1 = S_spare1;
othersignals.S_spare2 = S_spare2;
othersignals.S_spare3 = S_spare3;
othersignals.S_ref = S_ref;
clear S_spare* S_ref
vv = fieldnames(othersignals);

Secg = nan(size(signals,1),12,size(signals,2));
for i = 1:length(vv);
    Secg(:,i,:) = othersignals.(vv{i});
end
% =
X = Secg./(repmat(range(Secg),[size(Secg,1) 1 1])) + repmat([0:-1:-(size(Secg,2)-1)],[size(Secg,1) 1 size(Secg,3)]);

GUI_Select_Windows(X);
% templates
tt = Template.Windows{1};tt(tt<1|tt>size(Secg,1))=[]; %% Template for beat
ttTW = Template.Windows{2};ttTW(ttTW<1|ttTW>size(Secg,1))=[]; %% Template for T wave
tt_noiseTP = Template.Windows{3};tt_noiseTP(tt_noiseTP<1|tt_noiseTP>size(Secg,1))=[]; %% T Template for noise
tt_sigQRS = Template.Windows{4};tt_sigQRS(tt_sigQRS<1|tt_sigQRS>size(Secg,1))=[]; %% T Template for signal (QRS)
clear Secg


% Align
do_alignment = 1;
if do_alignment
    
    D = nan(size(signals,2),length(vv));
    h = waitbar(0,'Measuring delay ...');
    for iv = 1:length(vv);
        
        x = othersignals.(vv{iv})(tt_sigQRS,:); % use QRS
        xref =  nanmedian(x,2);
        
        
        for i = 1:size(x,2);
            [c,tau] = xcorr(x(:,i),xref);
            d = -tau(c==max(c));
            D(i,iv) = d;
            
        end
        waitbar(i/length(vv),h);
    end
    Dref = round(nanmedian(D,2));
    close(h)
    
    %     othersignals_pre_a = othersignals;
    h = waitbar(0,'Aligning ECG ...');
    for iv = 1:length(vv);
        x = othersignals.(vv{iv});
        xref =  nanmedian(x,2);
        
        for i = 1:length(Dref);
            d = Dref(i);
            if d<0
                x(1:end+d,i) = [othersignals.(vv{iv})(-d+1:end,i)];
                x(end+d+1:end,i) = repmat(othersignals.(vv{iv})(end,i),[size(othersignals.(vv{iv})(end+d+1:end,i),1) 1]);
            elseif d>0
                x(d:end,i) = [othersignals.(vv{iv})(1:end-d+1,i)];
                x(1:d-1,i) = repmat(othersignals.(vv{iv})(1,i),[size(othersignals.(vv{iv})(1:d-1,i),1) 1]);
                
            end
        end
        othersignals.(vv{iv}) = x;
        clear x
        waitbar(i/length(vv),h);
    end
    %         signals_pre_a = signals;
    close(h)
    
    xs = signals;
    for i = 1:length(Dref);
        d = Dref(i);
        if d<0
            xs(1:end+d,i) = [signals(-d+1:end,i)];
            xs(end+d+1:end,i) = repmat(signals(end,i),[size(signals(end+d+1:end,i),1) 1]);
        elseif d>0
            xs(d:end,i) = [signals(1:end-d+1,i)];
            xs(1:d-1,i) = repmat(signals(1,i),[size(signals(1:d-1,i),1) 1]);
            
        end
    end
    signals = xs;
    clear xs
end

% correlation after alignment
Secg = nan(size(signals,1),12,size(signals,2));
for i = 1:length(vv);
    Secg(:,i,:) = othersignals.(vv{i});
end
Sm = nanmedian(Secg,3);
answer=inputdlg({'Threshold for Beat';'Threshold for QRS';'Threshold for TW'},'CC',1,{'0.8','0.9','0.8'});

CorrLim_QRS = str2double(answer{2});
CorrLim_beat = str2double(answer{1}); % becasue very long (3 beats)
CorrLim_TW = str2double(answer{3});
Corr_beat = nan(size(Secg,3),length(vv));
Corr_TW = nan(size(Secg,3),length(vv));
Corr_QRS = nan(size(Secg,3),length(vv));
h = waitbar(0,'Correlation ...');
for i = 1:length(vv)
    
    c = corrcoef([Sm(tt,i) squeeze(Secg(tt,i,:))]);
    cTW = corrcoef([Sm(ttTW,i) squeeze(Secg(ttTW,i,:))]);
    cQRS = corrcoef([Sm(tt_sigQRS,i) squeeze(Secg(tt_sigQRS,i,:))]);
    
    Corr_beat(:,i) = c(2:end,1);
    Corr_TW(:,i) = cTW(2:end,1);
    Corr_QRS(:,i) = cQRS(2:end,1);
    waitbar(i/length(vv),h);
    
end
close(h);
iibeat_ok = nanmedian(Corr_beat,2)>CorrLim_beat & nanmedian(Corr_TW,2)>CorrLim_TW & nanmedian(Corr_QRS,2)>CorrLim_QRS;
disp(['CorrLim. Beat=',num2str(CorrLim_beat),';  QRS=',num2str(CorrLim_QRS),'; TW=',num2str(CorrLim_TW)]);

disp(['Good Points (Beat) = ',num2str(sum(nanmedian(Corr_beat,2)>CorrLim_beat)),' (',num2str(100*nanmean(nanmedian(Corr_beat,2)>CorrLim_beat),3),'%)']);
disp(['Good Points (TW) = ',num2str(sum(nanmedian(Corr_TW,2)>CorrLim_TW)),' (',num2str(100*nanmean(nanmedian(Corr_TW,2)>CorrLim_TW),3),'%)']);
disp(['Good Points (QRS) = ',num2str(sum(nanmedian(Corr_QRS,2)>CorrLim_QRS)),' (',num2str(100*nanmean(nanmedian(Corr_QRS,2)>CorrLim_QRS),3),'%)']);
disp(' ------------------ ')
disp(['Good Points (Beat + TW + QRS) = ',num2str(sum(iibeat_ok)),' (',num2str(100*nanmean(iibeat_ok),3),'%)']);

signals = signals(tt,iibeat_ok);
for iv = 1:length(vv);
    othersignals.(vv{iv}) = othersignals.(vv{iv})(tt,iibeat_ok);
end
ParamSig.Label = ParamSig.Label(iibeat_ok);
xyz = xyz(:,iibeat_ok);
xyz_surfP = xyz_surfP(:,iibeat_ok);
ParamAnalysis.Removed_corr_before_alignment = sum(iibeat_ok==0);

% -
Template_ORI = Template;
for i = 1:4
    if ~isempty(Template.Windows{i})
        Template.Windows{i} = Template_ORI.Windows{i}-Template_ORI.Windows{1}(1)+1; %% Template for beat
        Template.Windows{i}(Template.Windows{i}<1)=[];
    end
end





geo.Cmesh.triangles = ParamOut.MESH.Faces;
geo.Cmesh.vertices = ParamOut.MESH.vertices;
geo.Cmesh.vertices_norm = ParamOut.MESH.Normals_number;
geo.Study_name = ParamOut.Study_name;
geo.xyz = xyz.';
geo.xyz_surfP = xyz_surfP.';
%

clear s Sm
% delete duplicates
h = waitbar(0,'Correlation ...');
dd = nan(size(signals,2));
for i = 1:size(signals,2);
    dd(i,:) = sum(abs(signals(1:10:end,i)*ones(1,size(signals,2))-signals(1:10:end,:)));
    dd(i,i)= nan;
    waitbar(i/size(signals,2),h);
end
iirep = sum(triu(dd==0))>0;
close(h)

%     c = corrcoef([S_ref_m S_ref]);m90c = mean(c(2:end,1)>0.90);
%     ca = corrcoef([S_ref_m S_ref_a]);m90ca = mean(ca(2:end,1)>0.90);
%     disp([num2str(m90c*100,3),'% of signals with cc>0.90 wrt median template']);
%     disp([num2str(m90ca*100,3),'% of signals with cc>0.90 wrt median template after alignment']);

% Initialization (if d=0 do not change)


% Reject nans
iiko = mean(isnan(signals))>0.5 | iirep;
disp(['Elininating ',num2str(sum(iiko)),' channels'])
ParamSig.Label(iiko) = [];
signals(:,iiko) = [];
geo.xyz(iiko,:) = [];
geo.xyz_surfP(iiko,:) = [];
othersignals.S_ref(:,iiko) = [];
othersignals.S_spare1(:,iiko) = [];
othersignals.S_spare2(:,iiko) = [];
othersignals.S_spare3(:,iiko) = [];
ParamAnalysis.Removed_nan = sum(iiko);
ParamAnalysis.N_final = size(signals,2);
x = nan(size(signals,1),12,size(signals,2));
for i = 1:length(vv);
    x(:,i,:) = othersignals.(vv{i});
end
S_ECG = nanmedian(x,3);
clear x
%%
Correlation_param.CorrLim_QRS =CorrLim_QRS;
Correlation_param.CorrLim_beat = CorrLim_beat; % becasue very long (3 beats)
Correlation_param.CorrLim_TW = CorrLim_TW;

Correlation_param.Corr_beat = Corr_beat;
Correlation_param.Corr_TW = Corr_TW;
Correlation_param.Corr_QRS = Corr_QRS;
Correlation_param.iibeat_ok = iibeat_ok;
clear Corr_* CorrLim_*
%
ii = find(PATHNAME=='\');
file_save_name = [PATHNAME,PATHNAME(ii(end-1)+1:end-1),'_signals'];
[namesave,pathsave] = uiputfile('*.mat','Save As',file_save_name);
file_save_name = [pathsave,namesave];
disp(['Saving as: ',file_save_name])
save([file_save_name(1:end-4),'_original'],'S','Secg','iibeat_ok','Correlation_param','Template*','S_ECG','othersignals','ParamAnalysis','Correlation_param');
clear Secg S othersignals c*
%%
disp('... SNR')
Stemp = signals;
P = abs(fft(Stemp)).^2;
f = [1:size(P,1)]/size(P,1)*1000;
iisig = f>5 & f<=40;
iinoise = f>40 & f<=100;
SNR_spect = nanmean(P(iisig,:))./nanmean(P(iinoise,:));
SNR = 10.*log10(SNR_spect);
clear Stemp SNR_spect

S_sigQRS = signals(Template.Windows{4},:);
S_sigTW = signals(Template.Windows{2},:);
S_noise = signals(Template.Windows{3},:);
SNR_int_QRS = 10.*log10(var(detrend(S_sigQRS))./var(detrend(S_noise)));
SNR_int_TW = 10.*log10(var(detrend(S_sigTW))./var(detrend(S_noise)));
clear S_noise S_sigTW S_sigQRS
%
do_pacing = 1;
fs = ParamSig.frequency;
if do_pacing
    sensthr = 0.5;
    do_control = 0;
    minCL = 120;
    if exist('Scs','var')
        [spikes] = find_pacing_spikes_mo(Scs(:,10),fs,sensthr,do_control,minCL);
    else
        [spikes] = find_pacing_spikes_mo(signals(:,10),fs,sensthr,do_control,minCL/1000*ParamSig.frequency);
    end
    spikes = spikes/ParamSig.frequency*1000;
else
    [spikes,iiSD] = spikes_detection_SR(signals,ParamSig);
end

% Delete last spike
if exist('Scs','var')
    L= fix(size(signals,2)/50);
    Spikes_GUI_standalone(Scs(:,[1:L:end]),ParamSig.frequency,spikes,[],0);
else
    L= fix(size(signals,2)/50);
    Spikes_GUI_standalone(signals(:,[1:L:end]),ParamSig.frequency,spikes,ParamSig.Label([1:L:end]),0);
end
spikes = DataFromGui.spikes;
spikes_samp = round(spikes/1000*ParamSig.frequency);
Twidth_replace_spk = 30;
Thresh_replace_spk = 5;
[S_nosp,~] = replace_spikes_mo(signals,spikes_samp,round(Twidth_replace_spk/1000*fs),Thresh_replace_spk);

% Sig_corr
disp('... sig corr')
sp = spikes_samp;sp((size(signals,1)-sp)<400/1000*ParamSig.frequency)=[];
if length(sp)>1
    T = [round(10/1000*ParamSig.frequency) : round(nanmedian(diff(sp))/1000*ParamSig.frequency)];
else
    T = [round(10/1000*ParamSig.frequency) : round(0.45*ParamSig.frequency)];
end
[sig_corr] = sig_corr_mo_all_beats(signals,sp,T);

%% Markers (AT in non-filtered - RT in filtered)
disp('... filtering')
xx = max(Template.Windows{4}/ParamSig.frequency*1000)-spikes;xx(xx<0)=[];
ParamIn.DTmax = min(xx);

xx = min(Template.Windows{2}/ParamSig.frequency*1000)-spikes;xx(xx<0)=[];
ParamIn.min_RT = min(xx);

xx = max(Template.Windows{2}/ParamSig.frequency*1000)-spikes;xx(xx<0)=[];
ParamIn.max_RT = min(xx);
ParamIn.do_control = 1;
ParamIn.frequency = fs;

BW = 25; % only low pass (signal too short)
signals_proc = butterworthfilter(S_nosp,fs,BW,'low');
ParamIn.BW = [0 BW];

BW = 100; % only low pass (signal too short)
signals_proc_AT = butterworthfilter(S_nosp,fs,BW,'low');
ParamIn.BW_AT = [0 BW];
Answer = questdlg('Do you want to use a Notch filter?');
if isequal(Answer,'Yes');
    wo = 50/(ParamSig.frequency/2);
    bw = wo/35;
    [b,a] = iirnotch(wo,bw);
    signals_proc_AT = filtfilt(b,a,signals_proc_AT);
    ParamIn.Notch_AT = true;
else
    ParamIn.Notch_AT = false;
end
ParamIn.DeltaT = 15;

Markers_TW = DTRT_mo_gui_v4(signals_proc,spikes,ParamIn);
Markers_QRS = DTRT_mo_gui_v4(signals_proc_AT,spikes,ParamIn);

params_from_TW = {'tTpeak','tPwave_int','tPw','tTend','tTend_defl','rt_Wyatt','rt_Alternative','iiTwpos','rt_up','isot','rt_down','ATw','ATw2'};
params_from_QRS = {'tSw','tRw','tQw','tQRSon','tQRSoff','QRSw_fw90','QRSw','dt','tdtM','tdtm','RatioI','ParamIn','m_function','Legend',...,
    'Tpeak_amp','Rw_amp','Sw_amp','Qw_amp','QRS_amp','QRS_area'};

for i = 1:length(params_from_QRS)
    Markers.(params_from_QRS{i}) = Markers_QRS.(params_from_QRS{i});
end
for i = 1:length(params_from_TW)
    Markers.(params_from_TW{i}) = Markers_TW.(params_from_TW{i});
end


%=
% filename_save = [PATHNAME,FILENAME(1:end-8)];
save([file_save_name],'signals','ParamSig','SNR*','spikes','Markers','signals_proc*','geo','sig_corr');


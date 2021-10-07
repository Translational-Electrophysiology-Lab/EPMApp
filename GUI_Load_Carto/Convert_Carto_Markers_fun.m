function [] = Convert_Carto_Markers_fun(signals,ParamSig,Template);

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
iibeat_okf = find(iibeat_ok);
if do_pacing
    sensthr = 0.5;
    do_control = 0;
    minCL = 120;
    if exist('Scs','var')
        [spikes] = find_pacing_spikes_mo(Scs(:,iibeat_okf(1)),fs,sensthr,do_control,minCL);
    else
        [spikes] = find_pacing_spikes_mo(signals(:,iibeat_okf(1)),fs,sensthr,do_control,minCL);
    end
    
else
    [spikes,iiSD] = spikes_detection_SR(signals(:,iibeat_okf),ParamSig);
end

% Delete last spike
if exist('Scs','var')
    L= fix(sum(iibeat_ok)/50);
    Spikes_GUI_standalone(Scs(:,[iibeat_okf(1:L:end)]),ParamSig.frequency,spikes,[],0);
    clear iibeat_okf
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
T = [round(10/1000*ParamSig.frequency) : round(nanmedian(diff(sp))/1000*ParamSig.frequency)];
[sig_corr] = sig_corr_mo_all_beats(signals,sp,T);

%% Markers (AT in non-filtered - RT in filtered)
disp('... filtering')
xx = max(Template.Windows{4})-spikes;xx(xx<0)=[];
ParamIn.DTmax = min(xx);

xx = min(Template.Windows{2})-spikes;xx(xx<0)=[];
ParamIn.min_RT = max(min(xx),200);

xx = max(Template.Windows{2})-spikes;xx(xx<0)=[];
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

disp(['Saving: ',filename_save])
save([filename_save],'signals','ParamSig','SNR*','spikes','Markers','signals_proc*','geo','sig_corr');

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


% Correlation to eliminate signals even before alignment
keyboard
figure,plot(S_ref);
% samples_corr = [910:1850];
samples_corr = [1:size(S_ref,1)];

cc_th = 0.75;
c  = corrcoef([nanmedian(S_ref(samples_corr,:),2) S_ref(samples_corr,:)]);
iiko = c(1,2:end)<cc_th;
%
ParamSig.Label(iiko) = [];
signals(:,iiko) = [];
xyz(:,iiko) = [];
xyz_surfP(:,iiko) = [];
S_ref(:,iiko) = [];
S_spare1(:,iiko) = [];
S_spare2(:,iiko) = [];
S_spare3(:,iiko) = [];
ParamAnalysis.Removed_corr_before_alignment = sum(iiko);



% Align
do_alignment = 1;
if do_alignment
    S_ref_m = nanmedian(S_ref,2);
    S_ref_a = S_ref; % if D=0 do not change
    D = nan(1,size(S_ref_a,2));
    for i = 1:size(S_ref_a,2);
        [c,tau] = xcorr(S_ref(:,i),S_ref_m);
        d = -tau(c==max(c));
        D(i) = d;
        if d<0
            S_ref_a(1:end+d,i) = [S_ref(-d+1:end,i)];
            S_ref_a(end+d+1:end,i) = repmat(S_ref(end,i),[size(S_ref(end+d+1:end,i),1) 1]);
        elseif d>0
            S_ref_a(d:end,i) = [S_ref(1:end-d+1,i)];
            S_ref_a(1:d-1,i) = repmat(S_ref(1,i),[size(S_ref(1:d-1,i),1) 1]);
        else
        end
    end
    figure,
%     plot(S_ref(:,i)),hold on,plot(S_ref_m,'r'),plot(S_ref_a(:,i),'k')
    
    c = corrcoef([S_ref_m S_ref]);m90c = mean(c(2:end,1)>0.90);
    ca = corrcoef([S_ref_m S_ref_a]);m90ca = mean(ca(2:end,1)>0.90);
    disp([num2str(m90c*100,3),'% of signals with cc>0.90 wrt median template']);
    disp([num2str(m90ca*100,3),'% of signals with cc>0.90 wrt median template after alignment']);
    
    % Initialization (if d=0 do not change)
    signals_A = signals;
    S_spare1_A = S_spare1;
    S_spare2_A = S_spare2;
    S_spare3_A = S_spare3;
    
    for i = 1:size(S_ref_a,2);
        d = D(i);
        if d<0
            signals_A(1:end+d,i) = [signals(-d+1:end,i)];
            signals_A(end+d+1:end,i) = repmat(signals(end,i),[size(signals(end+d+1:end,i),1) 1]);
            
            S_spare1_A(1:end+d,i) = [S_spare1(-d+1:end,i)];
            S_spare1_A(end+d+1:end,i) = repmat(S_spare1(end,i),[size(S_spare1(end+d+1:end,i),1) 1]);
            
            S_spare2_A(1:end+d,i) = [S_spare2(-d+1:end,i)];
            S_spare2_A(end+d+1:end,i) = repmat(S_spare2(end,i),[size(S_spare2(end+d+1:end,i),1) 1]);
            
            S_spare3_A(1:end+d,i) = [S_spare3(-d+1:end,i)];
            S_spare3_A(end+d+1:end,i) = repmat(S_spare3(end,i),[size(S_spare3(end+d+1:end,i),1) 1]);
            
        elseif d>0
            signals_A(d:end,i) = [signals(1:end-d+1,i)];
            signals_A(1:d-1,i) = repmat(signals(1,i),[size(signals(1:d-1,i),1) 1]);
            
            S_spare1_A(d:end,i) = [S_spare1(1:end-d+1,i)];
            S_spare1_A(1:d-1,i) = repmat(S_spare1(1,i),[size(S_spare1(1:d-1,i),1) 1]);
            
            S_spare2_A(d:end,i) = [S_spare2(1:end-d+1,i)];
            S_spare2_A(1:d-1,i) = repmat(S_spare2(1,i),[size(S_spare2(1:d-1,i),1) 1]);
            
            S_spare3_A(d:end,i) = [S_spare3(1:end-d+1,i)];
            S_spare3_A(1:d-1,i) = repmat(S_spare3(1,i),[size(S_spare3(1:d-1,i),1) 1]);
            
        else
        end
    end
end

geo.MESH = ParamOut.MESH;
geo.Study_name = ParamOut.Study_name;
geo.xyz = xyz;
geo.xyz_surfP = xyz_surfP;
%
signals = signals_A;
othersignals.S_spare1 = S_spare1_A;
othersignals.S_spare2 = S_spare2_A;
othersignals.S_spare3 = S_spare3_A;
othersignals.S_ref = S_ref_a;
clear S_*
% Refine based on morphological correlations and nans
cc_th2 = 0.75;
vv = fieldnames(othersignals);
cc = nan(4,size(signals,2));
for iv = 1:length(vv)
   s = othersignals.(vv{iv});
   c = corrcoef([nanmedian(s,2) s]);
   cc(iv,:) = c(1,2:end);
end
iiko = sum(cc<cc_th2)>=1|mean(isnan(signals))>0.5;
ParamSig.Label(iiko) = [];
% ParamSig.
signals(:,iiko) = [];
geo.xyz(:,iiko) = [];
geo.xyz_surfP(:,iiko) = [];
othersignals.S_ref(:,iiko) = [];
othersignals.S_spare1(:,iiko) = [];
othersignals.S_spare2(:,iiko) = [];
othersignals.S_spare3(:,iiko) = [];
ParamAnalysis.Removed_corr_after_alignment = sum(iiko);
ParamAnalysis.N_final = size(signals,2);

%%
ii = find(PATHNAME=='\');
file_save_name = [PATHNAME,PATHNAME(ii(end-1)+1:end-1),'_signals'];
[namesave,pathsave] = uiputfile('*.mat','Save As',file_save_name);
file_save_name = [pathsave,namesave];
disp(['Saving as: ',file_save_name])
save(file_save_name,'geo','signals','ParamSig','othersignals','ParamAnalysis')

%%


function [signals,ParamSig,ParamOut] = Convert_Carto_fun(filename,do_pacing,ParamWin,do_markers)
if nargin<4
    do_markers = 0;
end
if nargin<3
    ParamWin.do_window = 1;
    ParamWin.cc_QRS = 0.9;
    ParamWin.cc_sig = 0.8;
    ParamWin.cc_TW = 0.8;
    
end
if nargin<2
    do_pacing = 0;
end


addpath E:\UCL\Scripts_all\Scripts_mo\GUI_egm\GUIs\
addpath E:\UCL\Scripts_all\Scripts_mo\GUI_egm\GUI_egm_mFiles\
addpath E:\UCL\Scripts_all\Scripts_mo\VT_RVI_Carto

ECG_labels = {'I(110)','II(111)','III(112)','aVL(171)','aVR(172)','aVF(173)','V1(22)','V2(23)','V3(24)','V4(25)','V5(26)','V6(27)'};
ParamSig.frequency = 1000;

%% Load mesh
[Cmesh] = LoadCartoMesh_fun([filename,'.mesh'],0);
geo.Cmesh = Cmesh.MESH;
geo.name = Cmesh.name;
%%
[CAR,Table_car] = LoadCartoCar_fun([filename,'_car.txt'],0);
geo.xyz = CAR.xyz;

catheter_IDs = [str2double(Table_car(:,20))];
catheter_type = unique(catheter_IDs);



%%
ii = find(filename=='\');
map_name = filename(ii(end)+1:end);
PATHNAME = filename(1:ii(end));
% SAVING files to analyse
dir_save = [PATHNAME,'MAT\'];
if ~exist(dir_save,'dir')
    mkdir(dir_save);
end
filename_save = [dir_save,map_name,'_',datestr(now,'yyyy_mm_dd')];

fileIDw = fopen([filename_save,'_EXPORT_LOG.txt'],'w');

for i = 1:length(catheter_type);
    fprintf(fileIDw,['n=',num2str(sum(catheter_IDs==catheter_type(i))),' points with catheter type = ',num2str(catheter_type(i)),'\r\n']);
end

cath_name = {'NAVISTAR_CONNECTOR','MAGNETIC_20_POLE_A_CONNECTOR','MAGNETIC_20_POLE_B_CONNECTOR','CS_CONNECTOR'};
% Get type of catheter
Xnavistar = dir([filename,'*NAVISTAR_CONNECTOR_Eleclectrode_Positions_OnAnnotation.txt']);
X20A = dir([filename,'*MAGNETIC_20_POLE_A_CONNECTOR_Eleclectrode_Positions_OnAnnotation.txt']);
X20B = dir([filename,'*MAGNETIC_20_POLE_B_CONNECTOR_Eleclectrode_Positions_OnAnnotation.txt']);
XCS = dir([filename,'*CS_CONNECTOR_Eleclectrode_Positions_OnAnnotation.txt']);

fprintf(fileIDw,'%s\r\n',[' Files - Navistar = ',num2str(length(Xnavistar))]);
fprintf(fileIDw,'%s\r\n',[' Files - 20A = ',num2str(length(X20A))]);
fprintf(fileIDw,'%s\r\n',[' Files - 20B = ',num2str(length(X20B))]);
fprintf(fileIDw,'%s\r\n',[' Files - CS = ',num2str(length(XCS))]);


ii = [isempty(Xnavistar) isempty(X20A) isempty(X20B) isempty(XCS)];
cath_name(ii) = [];

file_type = 'Eleclectrode_Positions_OnAnnotation';
% do_CS = 1; % To get precise localization of pacing spikes when pacing is on

Secg = nan(2500,12,length(CAR.Index_point));
S = cell(1,length(CAR.Index_point));
S_xyz = S;
Label_temp = cell(1,length(CAR.Index_point));
ParamSig.N_Exported_per_point = zeros(1,length(CAR.Index_point));
h = waitbar(0,[map_name,' Loading data ...']);


% check if CS signals are there
Xecg = dir([PATHNAME,map_name,'*export.txt']);
ii = randi(length(Xecg),1,10);
jjcs = zeros(1,length(ii));
for i = 1:length(ii)
    [~,ll] = LoadWaveformsFromCarto([PATHNAME,Xecg(ii(i)).name]);
    jjcs(i) = sum(contains(ll,'CS1'));
end
clear i ii ll
if sum(jjcs)>0
    do_CS = 1;
    Scs = nan(2500,length(CAR.Index_point));
else
    do_CS = 0;
    Scs = [];
end



% CathType20A = questdlg('Is 20A connected to a Decapolar catheter?','CATHETER','DECA','PENTA','PENTA');

cath_thype_all = cell(1,length(CAR.Index_point));
for ip = 1:length(CAR.Index_point);
    if rem(ip,fix(length(CAR.Index_point)/10))==0
        waitbar(ip/length(CAR.Index_point),h);
    end
    
    for icath = 1:length(cath_name)
        eln = [];
        loadname_pos = [PATHNAME,map_name,'_P',num2str(CAR.Index_point(ip)),'_',cath_name{icath},'_',file_type,'.txt'];
        % ==
        fid=fopen(loadname_pos);
        if fid==-1 % error with file or file does not exist
            continue
        end
        fgetl(fid);fgetl(fid);
        dataArray = textscan(fid,'%f %f %f %f %f','CollectOutput',1);
        fclose(fid);
        
        % ==
        if isequal(cath_name{icath},'MAGNETIC_20_POLE_A_CONNECTOR')|isequal(cath_name{icath},'MAGNETIC_20_POLE_B_CONNECTOR')|isequal(cath_name{icath},'CS_CONNECTOR') % Penta-array/Deca
            if size(dataArray{1},1)==22
                CathType20A = 'PENTA';
            elseif size(dataArray{1},1)==10
                CathType20A = 'DECA';
            else
                P = [dataArray{1}(:,3:end)];
                d = diff(sqrt(sum((P-P(1,:)).^2,2)));
                if d<2
                    CathType20A = 'PENTA';
                else
                    CathType20A = 'DECA';
                end
                
            end
            
            eln = find(sum(dataArray{1}(:,3:end) - ones(size(dataArray{1},1),1)*CAR.xyz(ip,:)==0,2)==3);
            eln = fix((eln-.01)/2)*2+[1:2]; % for each point, take the pair
            if ~isempty(eln)
                S_xyz{ip} =  dataArray{1}(eln,3:end);
                %                 if dataArray{1}(eln(1),3:end)~=CAR.xyz(ip,:);
                %                     error('Double check code here')
                %                 end
                if isequal(CathType20A,'PENTA')
                    eln = eln-2; % (The first 2 electrodes are sensors in Penta)
                end
            end
            
        else % Navistar
            eln = find(sum(dataArray{1}(:,3:end) - ones(size(dataArray{1},1),1)*CAR.xyz(ip,:)==0,2)==3); %
            eln = fix((eln-.01)/2)*2+[1:2]; % for each point, take the pair
            
            if isempty(eln)
                %                 load sensor and take the 4 poles
                loadname_sensor = [PATHNAME,map_name,'_P',num2str(CAR.Index_point(ip)),'_',cath_name{icath},'_Sensor_Positions_OnAnnotation.txt'];
                % ==
                fid=fopen(loadname_sensor);
                if fid==-1 % error with file
                    continue
                end
                fgetl(fid);fgetl(fid);
                dataArrayS = textscan(fid,'%f %f %f %f %f','CollectOutput',1);
                fclose(fid);
                if ~isempty(find(sum(dataArrayS{1}(:,3:end) - ones(size(dataArrayS{1},1),1)*CAR.xyz(ip,:)==0,2)==3))
                    eln = 1:4; % The
                end
            end
            S_xyz{ip} =  dataArray{1}(eln,3:end);
            
        end
        
        if ~isempty(eln) % stop if point has been identified
            break
        end
        clear dataArray loadname_pos
    end
    
    % =
    if isempty(eln)|fid==-1
        continue
    end
    load_sig = [PATHNAME,map_name,'_P',num2str(CAR.Index_point(ip)),'_ECG_Export.txt'];
    [signals,Labels,Info] = LoadWaveformsFromCarto(load_sig);
    % =
    
    if  isequal(cath_name{icath},'MAGNETIC_20_POLE_A_CONNECTOR')
        j = find(contains(Labels,'20A_1('));
        cath_thype_all{ip} = CathType20A;
    elseif isequal(cath_name{icath},'NAVISTAR_CONNECTOR')
        j = find(contains(Labels,'M1(1)'));
        cath_thype_all{ip} = 'ABL';
    elseif isequal(cath_name{icath},'MAGNETIC_20_POLE_B_CONNECTOR')
        j = find(contains(Labels,'20B_1('));
        cath_thype_all{ip} = CathType20A;
    elseif isequal(cath_name{icath},'CS_CONNECTOR')
        j = find(contains(Labels,'CS1('));
        cath_thype_all{ip} = CathType20A;
    else
        
        warning('Connector not found')
        continue
        
    end
    ic = j+eln-1;
    
    % =
    %     S(:,ip) = signals(:,ic);
    %     ParamSig.Label(ip) = Labels(ic);
    %     [~,iiECG] = intersect(Labels,ECG_labels,'stable');
    %     Secg(:,:,ip) = signals(:,iiECG);
    S{ip} = signals(:,ic);
    Label_temp{ip} = Labels(ic);
    [~,iiECG] = intersect(Labels,ECG_labels,'stable');
    Secg(:,:,ip) = signals(:,iiECG);
    
    %
    if do_CS
        ics = find(~cellfun(@isempty,strfind(Labels,'CS1(')));
        if isempty(ics)
            ics = find(~cellfun(@isempty,strfind(Labels,'CS1-')));
        end
        if ~isempty(ics)
            Scs(:,ip) = signals(:,ics);
        end
    end
    
    clear signals Labels
end
close(h);
Scs_ori = Scs;
Secg_ori = Secg;

cath_thype_all(cellfun(@isnumeric,cath_thype_all,'un',1))=[];

fprintf(fileIDw,'%s\r\n',['Points ABL CATH = ',num2str(sum(contains(cath_thype_all,'ABL')))]);
fprintf(fileIDw,'%s\r\n',['Points PENTARAY CATH = ',num2str(sum(contains(cath_thype_all,'PENTA')))]);
fprintf(fileIDw,'%s\r\n',['Points DECANAV CATH = ',num2str(sum(contains(cath_thype_all,'DECA')))]);



% Templates
X = Secg_ori./(repmat(range(Secg_ori),[size(Secg_ori,1) 1 1])) + repmat([0:-1:-11],[size(Secg_ori,1) 1 size(Secg_ori,3)]);

% Customise windows
if ParamWin.do_window
    GUI_Select_Windows_CARTO(X);
    cc = findobj(0,'name','GUI_Load_Carto');
    hmain = guidata(cc);
    Template = hmain.mat.Template;
    clear hmain
else
    Template.point = round(size(X,3));
    Template.Windows{1} = [1:2500];
    Template.Windows{2} = [2281:2450];
    Template.Windows{4} = [1901:2100];
    Template.Windows{3} = [1800:1900];
end

tt = Template.Windows{1};tt(tt<1|tt>size(Secg,1))=[]; %% Template for beat
ttTW = Template.Windows{2};ttTW(ttTW<1|ttTW>size(Secg,1))=[]; %% Template for T wave
tt_noiseTP = Template.Windows{3};tt_noiseTP(tt_noiseTP<1|tt_noiseTP>size(Secg,1))=[]; %% T Template for noise
tt_sigQRS = Template.Windows{4};tt_sigQRS(tt_sigQRS<1|tt_sigQRS>size(Secg,1))=[]; %% T Template for signal (QRS)


%%

if ~isempty(Scs)
    Scs = Scs_ori(tt,:);
end
%% Align based on spikes in CS only
do_align_spikes = do_pacing;
if do_align_spikes&~isempty(Scs)
    h2 = waitbar(0,'Aligning on spikes ...');
    sp = cell(1,size(Scs,2));
    for i = 1:size(Scs,2)
        h2 = waitbar(i/size(Scs,2),h2);
        try
            sp{i} = spikes_detection_SR(Scs(:,i),ParamSig);
            sp{i} = round(sp{i}/1000*ParamSig.frequency);
        catch ME
            disp('Spikes not found')
        end
        %         figure(1)
        %         hold off
        %         plot(Scs(:,i));
        %         x = round(sp{i}/1000*ParamSig.frequency);
        %         hold on
        %         plot(x,Scs(x,i),'or');
        %         title(num2str(i))
        %         pause
        
    end
    
    iiok = ~cellfun(@isempty,sp);
    s2 = nan(1,length(sp));
    s2(iiok) = cellfun(@(x) x(end),sp(iiok));
    sref = round(nanmedian(s2));
    
    fprintf(fileIDw,'%s\r\n',['Reference spike for alignment = ',num2str(sref), ' samples']);
    
    
    S_a = S;
    Secg_a = Secg;
    Scs_a = Scs;
    D = nan(1,length(S));
    Nspikes_min = 2; % Number of spikes in window [adapt this to each pt]
    for i = 1:length(S)
        if length(sp{i})>=Nspikes_min
            %             d = sp{i}(end)-sref;
            [~,ib] = min(abs(sp{i}-sref));
            d = sp{i}(ib)-sref; clear ib
            D(i) = d;
            if d>0
                Scs_a(:,i) = [Scs(d:end,i);Scs_a(1,i)*ones(d-1,1)];
                Secg_a(:,:,i) = [Secg_a(d:end,:,i);ones(d-1,1)*Secg_a(1,:,i)];
                S_a{i} = [S_a{i}(d:end,:);ones(d-1,1)*S_a{i}(1,:)];
            elseif d<0
                Scs_a(:,i) = [Scs_a(1,i)*ones(-d,1);Scs_a(1:end+d,i)];
                Secg_a(:,:,i) = [ones(-d,1)*Secg_a(1,:,i);Secg_a(1:end+d,:,i)];
                S_a{i} = [ones(-d,1)*S{i}(1,:);S{i}(1:end+d,:)];
            else
            end
        end
        
    end
    
    % correcting original ones
    S = S_a;
    Secg = Secg_a;
    Scs = Scs_a;
    clear S*_a
    close(h2)
end
% l = cellfun(@length,sp)==3;figure,plot(Scs_ori(tt,l),'k'),hold on,plot(Scs(:,l),'r')

%%
hf = figure;
ax(1)=subplot(221);
i = 2;
plot(squeeze(Secg(:,i,:)),'k'),hold on,plot(nanmedian(Secg(:,i,:),3),'-r'),hold on,plot(tt,Secg(tt,i, Template.point),'g')
title(ECG_labels{i})
% -
ax(2)=subplot(222);
i = 5;
plot(squeeze(Secg(:,i,:)),'k'),hold on,plot(nanmedian(Secg(:,i,:),3),'-r'),hold on,plot(tt,Secg(tt,i, Template.point),'g')
title(ECG_labels{i})
% -
ax(3)=subplot(223);
i = 8;
plot(squeeze(Secg(:,i,:)),'k'),hold on,plot(nanmedian(Secg(:,i,:),3),'-r'),hold on,plot(tt,Secg(tt,i, Template.point),'g')
title(ECG_labels{i})
% -
ax(4)=subplot(224);
i = 11;
plot(squeeze(Secg(:,i,:)),'k'),hold on,plot(nanmedian(Secg(:,i,:),3),'-r'),hold on,plot(tt,Secg(tt,i, Template.point),'g')
title(ECG_labels{i})
% -
axis(ax,'tight');
linkaxes(ax,'x')

if ParamWin.do_window==0
    do_median_template =1;
else
    ButtonName = questdlg('Which template sould be used?','Chose a template','Median (red)','Single beat (green)','Median (red)');
    if isequal(ButtonName,'Median (red)')
        do_median_template =1 ; % use median as template
    else
        do_median_template =0 ; % use single beat as template
    end
end


%% No alignment
% Sa = Secg;

%%
% answer=inputdlg({'Threshold for Beat';'Threshold for QRS';'Threshold for TW'},'CC',1,{'0.8','0.9','0.8'});

CorrLim_QRS = ParamWin.cc_QRS;
CorrLim_beat = ParamWin.cc_sig; % becasue very long (3 beats)
CorrLim_TW = ParamWin.cc_TW;

if do_median_template
    Sm = nanmedian(Secg,3);
else
    Sm = Secg(:,:,Template.point);
end
iileads_ok = true(size(ECG_labels));
if isfield(Template,'Lead_ko')
    [~,iileads] = intersect(ECG_labels,Template.Lead_ko);
    iileads_ok(iileads) = false;
    iileads = [];
end
Sm(:,~iileads_ok) = nan;
% =
Corr_beat = nan(size(Secg,3),12);
Corr_TW = nan(size(Secg,3),12);
Corr_QRS = nan(size(Secg,3),12);
h = waitbar(0,[map_name, ' Correlation and save ...']);
for i = 1:12
    if ~iileads_ok(i)
        continue
    end
    c = corrcoef([Sm(tt,i) squeeze(Secg(tt,i,:))]);
    cTW = corrcoef([Sm(ttTW,i) squeeze(Secg(ttTW,i,:))]);
    cQRS = corrcoef([Sm(tt_sigQRS,i) squeeze(Secg(tt_sigQRS,i,:))]);
    
    Corr_beat(:,i) = c(2:end,1);
    Corr_TW(:,i) = cTW(2:end,1);
    Corr_QRS(:,i) = cQRS(2:end,1);
    waitbar(i/14,h);
    
end
iibeat_ok = nanmedian(Corr_beat,2)>CorrLim_beat & nanmedian(Corr_TW,2)>CorrLim_TW & nanmedian(Corr_QRS,2)>CorrLim_QRS;
fprintf(fileIDw,['CorrLim. Beat=',num2str(CorrLim_beat),';  QRS=',num2str(CorrLim_QRS),'; TW=',num2str(CorrLim_TW),'\r\n']);

fprintf(fileIDw,'%s\r\n',['Good Points (Beat) = ',num2str(sum(nanmedian(Corr_beat,2)>CorrLim_beat)),' (',num2str(100*nanmean(nanmedian(Corr_beat,2)>CorrLim_beat),3),'%)']);
fprintf(fileIDw,'%s\r\n',['Good Points (TW) = ',num2str(sum(nanmedian(Corr_TW,2)>CorrLim_TW)),' (',num2str(100*nanmean(nanmedian(Corr_TW,2)>CorrLim_TW),3),'%)']);
fprintf(fileIDw,'%s\r\n',['Good Points (QRS) = ',num2str(sum(nanmedian(Corr_QRS,2)>CorrLim_QRS)),' (',num2str(100*nanmean(nanmedian(Corr_QRS,2)>CorrLim_QRS),3),'%)']);
fprintf(fileIDw,'%s\r\n',[' ------------------ ']);
fprintf(fileIDw,'%s\r\n',['Good Points (Beat + TW + QRS) = ',num2str(sum(iibeat_ok)),' (',num2str(100*nanmean(iibeat_ok),3),'%)']);

%

%% Discard bad signals and only keep the beat interval
S_ECG = nanmean(Secg(tt,1:12,iibeat_ok),3);
signals = [S{iibeat_ok}];
if isempty(signals)
    warning('No signals');
    close(h)
    return
end
signals = signals(tt,:);
% Secg_ori = Secg;
% Secg = Secg(tt,1:12,iibeat_ok);

xecg = Secg(tt,1:12,iibeat_ok);
S_ECG_beats = nan(size(signals,1),12,size(signals,2));
l = cellfun(@(x) size(x,2),S(iibeat_ok));
icsum = cumsum(l);icsum = [0 icsum];
for i = 1:sum(iibeat_ok)
    ii = icsum(i)+1 : icsum(i+1);
    S_ECG_beats(:,:,ii) = repmat(xecg(:,:,i),[1 1 length(ii)]);
end


S_xyz2 = cellfun(@(x) x.',S_xyz,'UniformOutput',0);
xyz = [S_xyz2{iibeat_ok}];

point_n = CAR.Index_point;
point_n(~iibeat_ok) = [];
L = Label_temp;
L(~iibeat_ok) = [];
ParamSig.Label = cell(1,size(signals,2));
ind = 0;
for i = 1:length(point_n);
    for j = 1:length(L{i});
        ind = ind+1;
        p = ['P',num2str(point_n(i)),'-',L{i}{j}];
        p([find(p=='('):find(p==')')])=[];
        ParamSig.Label{ind} = p;
    end
end

Template_ORI = Template;
for i = 1:4
    if ~isempty(Template.Windows{i})
        Template.Windows{i} = Template_ORI.Windows{i}-Template_ORI.Windows{1}(1)+1; %% Template for beat
        Template.Windows{i}(Template.Windows{i}<1)=[];
    end
end

%%
fprintf(fileIDw,'%s\r\n',['Exported ',num2str(size(signals,2)),' UNI Electrograms']);
fprintf(fileIDw,'%s\r\n',['Original points ',num2str(length(CAR.Index_point))]);

%
geo.xyz = xyz.';

do_eliminate_duplicate_position=1;
if do_eliminate_duplicate_position
    P = geo.xyz;
    Dmin = nan(1,size(P,1));
    iDmin = Dmin;
    Dzero = nan(size(P));
    for i = 1:size(P,1);
        d = sqrt(sum(( abs(ones(size(P,1),1)*P(i,:) - P).^2),2));
        d(i) = nan;
        [Dmin(i),iDmin(i)] = min(d);
        Dzero(:,i) = d==0;
    end
    DzeroU = triu(Dzero);
    iiko_pos = find(sum(DzeroU)>0);
    
    signals(:,iiko_pos) = [];
    geo.xyz(iiko_pos,:) = [];
    ParamSig.Label(iiko_pos) = [];
    S_ECG_beats(:,:,iiko_pos) = [];
end

% Amplitude analysis
do_amplitude_analysis=1;
if do_amplitude_analysis
    thr_amp = 4;
    iiko = range(signals)>thr_amp*nanmedian(range(signals));
    %     figure,plot(range(signals),'.-');hold on,plot([1 4000],[1 1]*nanmedian(range(signals))*thr_amp);
    signals(:,iiko) = [];
    geo.xyz(iiko,:) = [];
    ParamSig.Label(iiko) = [];
    S_ECG_beats(:,:,iiko) = [];
end


waitbar(13/14,h);
fprintf(fileIDw,'%s\r\n',['Saving ',filename_save]);
save(filename_save,'signals','ParamSig','map_name','geo','S_ECG_beats');

% - Correlations params
Correlation_param.CorrLim_QRS =CorrLim_QRS;
Correlation_param.CorrLim_beat = CorrLim_beat; % becasue very long (3 beats)
Correlation_param.CorrLim_TW = CorrLim_TW;

Correlation_param.Corr_beat = Corr_beat;
Correlation_param.Corr_TW = Corr_TW;
Correlation_param.Corr_QRS = Corr_QRS;
Correlation_param.iibeat_ok = iibeat_ok;

waitbar(13/14,h);
% SAVING files for tracking the analysis
filename_save_original = [filename_save,'_original'];
fprintf(fileIDw,'%s\r\n',['Saving: ',filename_save_original]);
save(filename_save_original,'S','Secg','Label_temp','iibeat_ok','Correlation_param','Template*','S_ECG','S_ECG_beats','CAR','ECG_labels');

% =
saveas(hf,filename_save,'fig')
close(hf)

clear ax
hf2 = figure;
ax(1) = subplot(221);
plot(squeeze(Secg_ori(:,2,:)),'k');
ax(2) = subplot(222);
plot(squeeze(Secg_ori(:,10,:)),'k');
ax(3) = subplot(223);
plot(squeeze(Secg(:,2,:)),'k');
hold on
plot(tt,squeeze(Secg(tt,2,iibeat_ok)),'--r');
ax(4) = subplot(224);
aa1=plot(squeeze(Secg(:,10,:)),'k');
hold on
aa2=plot(tt,squeeze(Secg(tt,10,iibeat_ok)),'--r');

set(ax,'xlim',[1 2500])
title(ax(1),['Original beats ',ECG_labels{2}])
title(ax(3),['Re-aligned beats ',ECG_labels{2}])
title(ax(2),['Original beats ',ECG_labels{10}])
title(ax(4),['Re-aligned beats ',ECG_labels{10}])
axis(ax,'tight')
ll = legend([aa1(1) aa2(1)],'ALL','OK');
set(ll,'location','eastoutside');
set(hf2,'units','centimeters','paperUnits','centimeters','position',[2 2 18 15],'paperposition',[0 0 18 15],'papersize',[ 18 15]);

set(ax(3),'position',[.08 .1 .355 .36])
set(ax(4),'position',[.5 .1 .355 .36])
set(ax(1),'position',[.08 .57 .355 .36])
set(ax(2),'position',[.5 .57 .355 .36])
set(ll,'position',[.87 .35 .12 .07])
set(ax,'fontsize',12)

xlabel(ax(3),'samples')
xlabel(ax(4),'samples')

% =
fprintf(fileIDw,'%s\r\n',['Saving: ',filename_save,'_Beat_Selection']);
saveas(hf2,[filename_save,'_Beat_Selection'],'fig')
close(hf2);
waitbar(1,h);
close(h)

clearvars -except signals Template Param* iibeat_ok filename* fileIDw do* geo Scs*

%% Markers
if do_markers
    disp('... SNR')
    Stemp = signals;
    P = abs(fft(Stemp)).^2;
    f = [1:size(P,1)]/size(P,1)*1000;
    iisig = f>5 & f<=40;
    iinoise = f>40 & f<=100;
    SNR_spect = nanmean(P(iisig,:))./nanmean(P(iinoise,:));
    SNR = 10.*log10(SNR_spect);
    clear Stemp SNR_spect
    
%     S_sigQRS = signals(Template.Windows{4},:);
%     S_sigTW = signals(Template.Windows{2},:);
%     clear  S_sigTW S_sigQRS
    %
%     do_pacing = 1;
    fs = ParamSig.frequency;
    iibeat_okf = find(iibeat_ok);
    do_control = 0;
    if do_pacing
        sensthr = 0.5;
        minCL = 120;
                   [spikes] = find_pacing_spikes_mo(signals(:,iibeat_okf(1)),fs,sensthr,do_control,minCL);
        
    else
            [spikes] = spikes_detection_SR_CinC(signals(:,iibeat_okf(1)),fs);
%         [spikes] = find_pacing_spikes_mo(signals(:,iibeat_okf(1)),fs,sensthr,do_control,minCL);
    end
    
    % Delete last spike
    if ~isempty(Scs)&&do_pacing
        L= fix(sum(iibeat_ok)/50);
        if L==0
            L=1;
        end
        Spikes_GUI_standalone_CARTO_check(Scs(:,[iibeat_okf(1:L:end)]),ParamSig.frequency,spikes,[],0);
        clear iibeat_okf
    else
        L= fix(size(signals,2)/50);
        Spikes_GUI_standalone_CARTO_check(signals(:,[1:L:end]),ParamSig.frequency,spikes,ParamSig.Label([1:L:end]),0);
    end
    
    cc = findobj(0,'name','GUI_Load_Carto');
    hmain = guidata(cc);
    spikes = hmain.mat.spikes;
    clear hmain
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
        T = [round(10/1000*ParamSig.frequency) : round(450/1000*ParamSig.frequency)];
    end
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
    
    ParamIn.Notch_AT = false;
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
    
    fprintf(fileIDw,'%s\r\n',['Saving Markers in ',filename_save]);
    disp(['Saving: ',filename_save])
    save([filename_save],'signals','ParamSig','SNR*','spikes','Markers','signals_proc*','geo','sig_corr');
    
    
end
fclose(fileIDw);

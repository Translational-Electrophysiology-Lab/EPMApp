clear

% dir_main = 'F:\Data_Barts_VT_ECGi_ENDO\ORI\';
% dir_save_main = 'F:\Data_Barts_VT_ECGi_ENDO\MAT\';

dir_main = uigetdir; 
dir_main = [dir_main,'\'];
% dir_main = 'E:\UCL\Data&People\Data_ECGI_Shocks\ECGI_Risk_Export\';
dir_save_main = [dir_main(1:end-1),'_MAT\'];
if ~exist(dir_save_main,'dir')
    mkdir(dir_save_main)
end

do_plot = 1;
%
X = dir([dir_main]);
X(strncmp({X.name},'.',1))=[];
X(cellfun(@(x) x==0,{X.isdir}))=[];
[~,is] = sort({X.date});
sbj_name_all_ori = {X(is(end:-1:1)).name};

Xdone = dir([dir_save_main]);
Xdone(strncmp({Xdone.name},'.',1))=[];
Xdone(cellfun(@(x) x==0,{Xdone.isdir}))=[];
sbj_name_done = {Xdone.name};

ii = listdlg('ListString',sbj_name_all_ori,'SelectionMode','multiple');
sbh_to_do = sbj_name_all_ori(ii);
[sbj_name_done_s,iidone] = intersect(sbh_to_do,sbj_name_done,'stable');

if ~isempty(iidone)
    Ans = questdlg(['These pts have already been converted: discard?',sbj_name_done_s]);
    if isequal(Ans,'Yes')
        sbh_to_do(iidone) = [];
    end
end

for isbj = 1 : length(sbh_to_do)
    
    %% Geometry
    dir_geometry = [dir_main,sbh_to_do{isbj},'\GEO\'];
    
    Xfilename = dir([dir_geometry,'*.Ventricles.Vertices.csv']);
    Xfilename(strncmp({Xfilename.name},'.',1))=[];
    
    jj = strfind(Xfilename(1).name,'Ventricles.Vertices.csv');
    filename = Xfilename(1).name(1:jj-1);
    filename_save = sbh_to_do{isbj};
    
    
    % == Load geometry [Ventricles]
    if exist([dir_geometry,filename,'Ventricles.Faces.csv'],'file')
        T = csvread([dir_geometry,filename,'Ventricles.Faces.csv']);
        geo.Mesh_ventricles = T+1;
        T = csvread([dir_geometry,filename,'Ventricles.Vertices.csv']);
        geo.Nodes_ventricles = T*10; % cm -> mm (double check)
        ParamSig.dir_geo = dir_geometry;
    else
        warning(['Missing Ventricles']);
    end
    
    % == Load geometry [Aorta]
    if exist([dir_geometry,filename,'Aorta_Landmark.Faces.csv'],'file')
        T = csvread([dir_geometry,filename,'Aorta_Landmark.Faces.csv']);
        geo.Mesh_Aorta = T+1;
        T = csvread([dir_geometry,filename,'Aorta_Landmark.Vertices.csv']);
        geo.Nodes_Aorta = T*10; % cm -> mm (double check)
    else
        warning(['Missing Aorta']);
    end
    
    % == Load geometry [Valves]
    if exist([dir_geometry,filename,'Valves_Landmark.Faces.csv'],'file')
        T = csvread([dir_geometry,filename,'Valves_Landmark.Faces.csv']);
        geo.Mesh_Valves = T+1;
        T = csvread([dir_geometry,filename,'Valves_Landmark.Vertices.csv']);
        geo.Nodes_Valves = T*10; % cm -> mm (double check)
    else
        warning(['Missing Valves']);
    end
    
    % == Load geometry [LAD]
    if exist([dir_geometry,filename,'LAD_Landmark.Faces.csv'],'file')
        T = csvread([dir_geometry,filename,'LAD_Landmark.Faces.csv']);
        geo.Mesh_LAD = T+1;
        T = csvread([dir_geometry,filename,'LAD_Landmark.Vertices.csv']);
        geo.Nodes_LAD = T*10; % cm -> mm (double check)
    else
        warning(['Missing LAD']);
    end
    
    % == Load geometry [Torso]
    if exist([dir_geometry,filename,'Torso.Faces.csv'],'file')
        T = csvread([dir_geometry,filename,'Torso.Faces.csv']);
        geo.Mesh_Torso = T+1;
        T = csvread([dir_geometry,filename,'Torso.Vertices.csv']);
        geo.Nodes_Torso = T*10; % cm -> mm (double check)
    else
        warning(['Missing Torso']);
    end
    
    % == Load geometry [Atria]
    if exist([dir_geometry,filename,'Atria_Landmark.Faces.csv'],'file')
        T = csvread([dir_geometry,filename,'Atria_Landmark.Faces.csv']);
        geo.Mesh_Atria = T+1;
        T = csvread([dir_geometry,filename,'Atria_Landmark.Vertices.csv']);
        geo.Nodes_Atria = T*10; % cm -> mm (double check)
        ParamSig.dir_geo = dir_geometry;
    else
        warning(['Missing Atria']);
    end
    
        % == Load geometry [LV lead]
    if exist([dir_geometry,filename,'lv lead_Landmark.Faces.csv'],'file')
        T = csvread([dir_geometry,filename,'lv lead_Landmark.Faces.csv']);
        geo.Mesh_LV_lead = T+1;
        T = csvread([dir_geometry,filename,'lv lead_Landmark.Vertices.csv']);
        geo.Nodes_LV_lead = T*10; % cm -> mm (double check)
        ParamSig.dir_geo = dir_geometry;
    else
        warning(['Missing LV Lead']);
    end
    
     % == Load geometry [LV lead]
    if exist([dir_geometry,filename,'rv lead_Landmark.Faces.csv'],'file')
        T = csvread([dir_geometry,filename,'rv lead_Landmark.Faces.csv']);
        geo.Mesh_RV_lead = T+1;
        T = csvread([dir_geometry,filename,'rv lead_Landmark.Vertices.csv']);
        geo.Nodes_RV_lead = T*10; % cm -> mm (double check)
        ParamSig.dir_geo = dir_geometry;
    else
        warning(['Missing RV Lead']);
    end

    
    
    dir_sig = [dir_main,sbh_to_do{isbj},'\MAP\'];
    Xsig = dir([dir_sig,'*_Potentials.txt']);
    Xsig(strncmp({Xsig.name},'.',1))=[];
    
    dir_save = [dir_save_main,sbh_to_do{isbj},'\'];
    if ~exist(dir_save,'dir')
        mkdir(dir_save);
    end
    
    if do_plot
        fn = fieldnames(geo);
        ii = find(contains(fn,'Nodes_')&~contains(fn,'Torso'));
        h= gobjects(1,length(ii));
        hf = figure;
        ll = lines;
        for ip = 1:length(ii)
            hold on
            xyz = geo.(fn{ii(ip)});
            tri = geo.(['Mesh_',fn{ii(ip)}(7:end)]);
            h(ip) = trisurf(tri,xyz(:,1),xyz(:,2),xyz(:,3));
            set(h(ip),'facecolor',ll(ip,:));
        end
        axis(gca,'image','off')
        set(h,'Edgecolor','none');
        camlight('head')
        title(sbh_to_do{isbj})
        
        saveas(hf,[dir_save,sbh_to_do{isbj},'_Anatomy'],'fig');
        close(hf)
    end
    
    
    for jsig = 1:length(Xsig)
        jj = strfind(Xsig(jsig).name,'_Potentials.txt');
        sig_name = [filename_save,'_',Xsig(jsig).name(1:jj-1)];
        
        % Load intracardiac
        fileID = fopen([Xsig(jsig).folder,'\',Xsig(jsig).name],'r');
        a = fgetl(fileID);
        N = sum(a==',')+1;
        fseek(fileID,0,-1);
        formatSpec = repmat('%f',[1 N]);
        dataArray = textscan(fileID, formatSpec, 'Delimiter', ',', 'EmptyValue' ,NaN, 'ReturnOnError',false,'collectoutput',1);
        fclose(fileID);
        signals = [dataArray{1}];
        clear dataArray
        % Load body surface
        fileID = fopen([Xsig(jsig).folder,'\',Xsig(jsig).name(1:jj-1),'_InputSignal.txt'],'r');
        a = fgetl(fileID);
        N = sum(a==',')+1;
        fseek(fileID,0,-1);
        formatSpec = repmat('%f',[1 N]);
        dataArray = textscan(fileID, formatSpec, 'Delimiter', ',', 'EmptyValue' ,NaN, 'ReturnOnError',false,'collectoutput',1);
        fclose(fileID);
        signals_body_surface = [dataArray{1}];
        clear dataArray
        
        
        ParamSig.frequency = 1000;
        ParamSig.geoname = 'CardioInsight';
        
        filename(filename==' '|filename=='.'|filename=='-') = '_';
        sig_name(sig_name==' '|sig_name=='.'|sig_name=='-') = '_';
        
        
        % Load AT
        fAT = [Xsig(jsig).folder,'\',Xsig(jsig).name(1:jj-1),'_ACtivation.txt'];
        if exist(fAT,'file')
            fileID = fopen(fAT,'r');
            a = fgetl(fileID);
            N = sum(a==',')+1;
            fseek(fileID,0,-1);
            formatSpec = repmat('%f',[1 N]);
            dataArray = textscan(fileID, formatSpec, 'Delimiter', ',', 'EmptyValue' ,NaN, 'ReturnOnError',false,'collectoutput',1);
            fclose(fileID);
            MarkersC.dt = [dataArray{1}-1]/ParamSig.frequency*1000;
            MarkersC.rt_Wyatt  = nan(size(MarkersC.dt ));
            MarkersC.rt_Alternative  = nan(size(MarkersC.dt ));
            
            spikes = 1/ParamSig.frequency*1000;
        end
        clear dataArray
    
        
        
        %% saving
        %         dir_save = [dir_save_main,filename(1:end-1),'\'];
        
        if ~isdir(dir_save)
            mkdir(dir_save);
        end
        filename_intra = [dir_save,sig_name,'_Intracardiac'];
        filename_bodysurf = [dir_save,sig_name,'_BodySurface'];
        
        disp(['Saving ',filename_intra]);
        if exist(fAT,'file')
            save(filename_intra,'signals','geo','ParamSig','spikes','MarkersC');
        else
            save(filename_intra,'signals','geo','ParamSig');
        end
        clear signals
        % =
        disp(['Saving ',filename_bodysurf]);
        signals = signals_body_surface;
        save(filename_bodysurf,'signals','geo','ParamSig');
        clear signals
        
    end
    clear geo ParamSig
end


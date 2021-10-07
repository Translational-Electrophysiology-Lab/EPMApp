function OUT = Load_case_Precision_fun_for_GUI(meshname)

if  isequal(meshname(end-3:end),'.zip') % zipped folder
    
    h = helpdlg('Please wait while unzipping')
    f0 = unzip(meshname,meshname(1:end-4));
    delete(h)
    ii = find(contains(f0,'DxLandmarkGeo.xml'));
    str = cell(1,length(ii));
    for j = 1:length(ii)
        p = fileparts(f0{ii(j)});
        ii0 = find(p==filesep);
        str{j} = p(ii0(end)+1:end);
    end
    indx = listdlg('ListString',str,'PromptString','Please seelect one map','SelectionMode','single','ListSize',[300 200]);
    PATHNAME = [fileparts(f0{ii(indx)}),filesep];
    [ParamOut] = load_DxLandmarkGeo_Precision(f0{ii(indx)},0);
    meshname = [meshname(1:end-4),'\'];
    file_name = fileparts(f0{ii(indx)});
    clear ii
    ii = find(file_name==filesep);
    file_name = file_name(ii(end)+1:end);
    clear ii
else
    error('Please compress the data into a zip file')
end



OUT.MESH = ParamOut;

disp(['... loading files in ',PATHNAME])
% -
X = dir([PATHNAME,'DxL*.csv*']);


OUT.file_name = file_name;
OUT.S = [];
OUT.labels = {};
OUT.xyz = []; OUT.xyz_surfP = [];
OUT.S_ref = [];OUT.S_spare1 = [];
OUT.S_spare2 = [];OUT.S_spare3 = [];
OUT.S_spare1_name = {};
OUT.S_spare2_name = {};
OUT.S_spare3_name = {};

OUT.Utilized = [];
OUT.Amp_neg = [];
OUT.Amp_p2p = [];
OUT.Time_LAT_ref = [];
OUT.Time_LAT_rov = [];
OUT.Time_End = [];
OUT.CFE_m = [];
OUT.CFE_sd = [];

OUT.Point_number = [];
OUT.CL = [];
OUT.Segment = {};
OUT.Map_name = {};
OUT.Rov_detect = {};
OUT.Ref_detect = {};

ii = find(PATHNAME=='\');

h = waitbar(0,'Converting data ...');
for ip = 1:length(X)
    waitbar(ip/length(X),h);
    filename = [PATHNAME,X(ip).name];
    
    if isequal(filename(end-2:end),'.gz')
        f0 = gunzip(filename);
        [s,p] = load_DxL_Precision_v2(f0{1});
        delete(f0{1})
    else
        [s,p] = load_DxL_Precision_v2(filename);
    end
    
    if ip == 1
        OUT.fs = p.fs;
    end
    OUT.S = [OUT.S s];
    OUT.labels = cat(2,OUT.labels,p.Label);
    OUT.Utilized = [OUT.Utilized,p.Utilized'];
    OUT.xyz = [OUT.xyz p.xyz'];
    OUT.xyz_surfP = [OUT.xyz_surfP p.xyz_surfP'];
    
    OUT.S_ref = [OUT.S_ref p.Other_signals.signals_ref];
    OUT.S_spare1 = [OUT.S_spare1 p.Other_signals.signals_spare1];
    OUT.S_spare2 = [OUT.S_spare2 p.Other_signals.signals_spare2];
    OUT.S_spare3 = [OUT.S_spare3 p.Other_signals.signals_spare3];
    
    OUT.S_spare1_name = cat(2,OUT.S_spare1_name,p.S1_name');
    OUT.S_spare2_name = cat(2,OUT.S_spare2_name,p.S2_name');
    OUT.S_spare3_name = cat(2,OUT.S_spare3_name,p.S3_name');
    
    OUT.Amp_neg = [OUT.Amp_neg,p.Amp_neg'];
    OUT.Amp_p2p = [OUT.Amp_p2p,p.Amp_p2p'];
    
    OUT.Time_LAT_ref = [OUT.Time_LAT_ref,p.Time_LAT_ref'];
    OUT.Time_LAT_rov = [OUT.Time_LAT_rov,p.Time_LAT_rov'];
    OUT.Time_End = [OUT.Time_End,p.Time_end'];
    
    OUT.CFE_m = [OUT.CFE_m,p.CFE_m'];
    OUT.CFE_sd = [OUT.CFE_sd,p.CFE_sd'];
    OUT.Point_number = [OUT.Point_number,p.Point_number'];
    OUT.CL = [OUT.CL,p.CL'];
    OUT.Segment = cat(1,OUT.Segment,p.Segment);
    OUT.Map_name = cat(1,OUT.Map_name,p.Map_name);
    OUT.Rov_detect = cat(2,OUT.Rov_detect,p.Rov_detect');
    OUT.Ref_detect = cat(2,OUT.Ref_detect,p.Ref_detect');
end
close(h);


OUT.Map_name = unique(OUT.Map_name);
OUT.Segment = unique(OUT.Segment);
OUT.N_points_all = size(OUT.S,2);
u = OUT.Utilized;
vv = fieldnames(OUT);
for i = 1:length(vv)
    if size(OUT.(vv{i}),2)==length(u)
        OUT.(vv{i})(:,u==0) = [];
    end
end

for i = length(f0):-1:1
    if isfolder(f0{i})
        rmdir(f0{i})
    else
        delete(f0{i})
    end
end
X = dir(meshname);X(~[X.isdir])=[];X(strncmp({X.name},'.',1))=[];
for i = 1:length(X)
    rmdir([meshname,X(i).name])
end
rmdir(meshname)
clear ii i

p = [meshname(1:end-1),'_MAT',filesep];
if ~exist(p,'dir')
    mkdir(p);
end

disp(['Saving all data in: ',[p,file_name]])
save([p,file_name],'OUT')

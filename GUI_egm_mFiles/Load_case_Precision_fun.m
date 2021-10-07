function OUT = Load_case_Precision_fun(meshname)

do_type = 0;
if isequal(meshname(end-2:end),'.gz') % gz files (all files zipped)
 do_type =1;   
    f0 = gunzip(meshname);
    [ParamOut] = load_DxLandmarkGeo_Precision(f0{1},0);
    delete(f0{1})
elseif isequal(meshname(end-3:end),'.zip') % zipped folder
    do_type=2;
    f0 = unzip(meshname,meshname(1:end-4));
    ii = find(contains(f0,'DxLandmarkGeo.xml'));
    [ParamOut] = load_DxLandmarkGeo_Precision(f0{ii},0);
    meshname = [meshname(1:end-4),'\'];

else
    do_type = 3;
    [ParamOut] = load_DxLandmarkGeo_Precision(meshname,0);
end
OUT.MESH = ParamOut;

ii = find(meshname==filesep);
PATHNAME = meshname(1:ii(end));
X = dir([PATHNAME,'DxL*.csv*']);



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
Tab.ID_map = PATHNAME(ii(end-1)+1:ii(end)-1);

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

% ParamAnalysis.N_ori = size(OUT.S,2);
% ParamAnalysis.Utilized_Ensite = sum(OUT.Utilized);
% ParamAnalysis.Removed_unused = sum(OUT.Utilized==0);
% ParamSig.frequency = OUT.fs;

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

if do_type==2
    for i = 1:length(f0)
        delete(f0{i})
    end
    rmdir(meshname)
    clear ii i
end
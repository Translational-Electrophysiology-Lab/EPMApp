function [signals,Param] = load_DxL_Precision_v2(filename);

% Example:
% filename = 'G:\Others\Velocity_Export\study_dwsG600573_2018_02_13_13_10_26\2018_03_16_17_57_16\DxL_14.csv';
% [signals,Param] = load_DxL_Precision(filename);
Param = [];
fid=fopen(filename);

%
str = 'Export from Study :';
a = '';
while ~contains(a,str)
    a = fgetl(fid);
    if a==-1;return;end
end
ii = find(a==':');
x = a(ii(end)+2:end);
Param.Study = x;

%
str = 'Export from Segment';
a = '';
while ~contains(a,str)
    a = fgetl(fid);
    if a==-1;return;end
end
ii = find(a==':');
x = a(ii(end)+2:end);
Param.Segment = x;

%
str = 'Export files Stored in Dir';
a = '';
while ~contains(a,str)
    a = fgetl(fid);
    if a==-1;return;end
end
ii = find(a==':');
x = a(ii(end)+2:end);
Param.Dir_save_ori = x;

%
str = 'Total number of data points (columns)';
a = '';
while ~contains(a,str)
    a = fgetl(fid);
    if a==-1;return;end
end
ii = find(a==':');
x = a(ii(end)+3:end);
Param.N = str2double(x);

%
str = 'This is file ';
a = '';
while ~contains(a,str)
    a = fgetl(fid);
    if a==-1;return;end
end
ii = find(a==',');
x = a(ii(end)+1:end);
Param.Map_name = x;
 
% Point number
str = 'pt number:';
a = '';
while ~contains(a,str)
    a = fgetl(fid);
    if a==-1;return;end
end
ii = find(a==',');
x = textscan(a(ii(1)+1:end),'%f','delimiter',',','collectoutput',1);
Param.Point_number = x{1};clear x

% ref trace name
str = 'ref trace:';
a = '';
while ~contains(a,str)
    a = fgetl(fid);
    if a==-1;return;end
end
ii = find(a==',');
if length(ii)>1
x = textscan(a(ii(1)+1:end),'%s','delimiter',',','collectoutput',1);
Param.Ref_trace_name = x{1};clear x
else
Param.Ref_trace_name = a(ii(1)+1:end);clear x    
end
% S1 trace name
str = 'S 1 trace:';
a = '';
while ~contains(a,str)
    a = fgetl(fid);
    if a==-1;return;end
end
ii = find(a==',');
x = textscan(a(ii(1)+1:end),'%s','delimiter',',','collectoutput',1);
Param.S1_name = x{1};clear x

% S2 trace name
str = 'S 2 trace:';
a = '';
while ~contains(a,str)
    a = fgetl(fid);
    if a==-1;return;end
end
ii = find(a==',');
x = textscan(a(ii(1)+1:end),'%s','delimiter',',','collectoutput',1);
Param.S2_name = x{1};clear x

% S3 trace name
str = 'S 3 trace:';
a = '';
while ~contains(a,str)
    a = fgetl(fid);
    if a==-1;return;end
end
ii = find(a==',');
x = textscan(a(ii(1)+1:end),'%s','delimiter',',','collectoutput',1);
Param.S3_name = x{1};clear x

% End Time
str = 'end time:';
a = '';
while ~contains(a,str)
    a = fgetl(fid);
    if a==-1;return;end
end
ii = find(a==',');
xend = textscan(a(ii(1)+1:end),'%f','delimiter',',','collectoutput',1);
Param.Time_end =  xend{1};
clear xend

% LAT
str = 'ref LAT:';
a = '';
while ~contains(a,str)
    a = fgetl(fid);
    if a==-1;return;end
end
ii = find(a==',');
xref = textscan(a(ii(1)+1:end),'%f','delimiter',',','collectoutput',1);

str = 'rov LAT:';
a = '';
while ~contains(a,str)
    a = fgetl(fid);
    if a==-1;return;end
end
ii = find(a==',');
xrov = textscan(a(ii(1)+1:end),'%f','delimiter',',','collectoutput',1);
Param.Time_LAT_rov =  xrov{1};
Param.Time_LAT_ref =  xref{1};
clear x*

% Amplitude, peak to peak
str = 'peak2peak:';
a = '';
while ~contains(a,str)
    a = fgetl(fid);
    if a==-1;return;end
end
ii = find(a==',');
x = textscan(a(ii(1)+1:end),'%f','delimiter',',','collectoutput',1);
Param.Amp_p2p = x{1};clear x

% Amplitude, peak neg
str = 'peak neg:';
a = '';
while ~contains(a,str)
    a = fgetl(fid);
    if a==-1;return;end
end
ii = find(a==',');
x = textscan(a(ii(1)+1:end),'%f','delimiter',',','collectoutput',1);
Param.Amp_neg = x{1};clear x

% Complex Fractionated Electrogran (CFE) mean (mv)
str = 'CFE mean:';
a = '';
while ~contains(a,str)
    a = fgetl(fid);
    if a==-1;return;end
end
ii = find(a==',');
x = textscan(a(ii(1)+1:end),'%f','delimiter',',','collectoutput',1);
Param.CFE_m = x{1};clear x

% Complex Fractionated Electrogran (CFE) sd (mv)
str = 'CFE stddev:';
a = '';
while ~contains(a,str)
    a = fgetl(fid);
    if a==-1;return;end
end
ii = find(a==',');
x = textscan(a(ii(1)+1:end),'%f','delimiter',',','collectoutput',1);
Param.CFE_sd = x{1};clear x

% XYZ roving
xyz = nan(Param.N,3);
str = 'roving x:';
a = '';
while ~contains(a,str)
    a = fgetl(fid);
    if a==-1;return;end
end
ii = find(a==',');
x = textscan(a(ii(1)+1:end),'%f','delimiter',',','collectoutput',1);
xyz(:,1) = x{1}(:);clear x
str = 'roving y:';
a = '';
while ~contains(a,str)
    a = fgetl(fid);
    if a==-1;return;end
end
ii = find(a==',');
x = textscan(a(ii(1)+1:end),'%f','delimiter',',');
xyz(:,2) = x{1}(:);clear x
str = 'roving z:';
a = '';
while ~contains(a,str)
    a = fgetl(fid);
    if a==-1;return;end
end
ii = find(a==',');
x = textscan(a(ii(1)+1:end),'%f','delimiter',',');
xyz(:,3) = x{1}(:);clear x
% =
xyz_surfP = nan(Param.N,3);
str = 'surfPt x:';
a = '';
while ~contains(a,str)
    a = fgetl(fid);
    if a==-1;return;end
end
ii = find(a==',');
x = textscan(a(ii(1)+1:end),'%f','delimiter',',');
xyz_surfP(:,1) = x{1}(:);clear x
str = 'surfPt y:';
a = '';
while ~contains(a,str)
    a = fgetl(fid);
    if a==-1;return;end
end
ii = find(a==',');
x = textscan(a(ii(1)+1:end),'%f','delimiter',',');
xyz_surfP(:,2) = x{1}(:);clear x
str = 'surfPt z:';
a = '';
while ~contains(a,str)
    a = fgetl(fid);
    if a==-1;return;end
end
ii = find(a==',');
x = textscan(a(ii(1)+1:end),'%f','delimiter',',');
xyz_surfP(:,3) = x{1}(:);clear x
% - 
str = 'cycle len:';
a = '';
while ~contains(a,str)
    a = fgetl(fid);
    if a==-1;return;end
end
ii = find(a==',');
x = textscan(a(ii(1)+1:end),'%f','delimiter',',');
Param.CL = x{1}(:);clear x
%
str = 'utilized:';
a = '';
while ~contains(a,str)
    a = fgetl(fid);
    if a==-1;return;end
end
ii = find(a==',');
x = textscan(a(ii(1)+1:end),'%f','delimiter',',');
Param.Utilized = x{1}(:);clear x
%
str = 'rov detect:';
a = '';
while ~contains(a,str)
    a = fgetl(fid);
    if a==-1;return;end
end
ii = find(a==',');
x = textscan(a(ii(1)+1:end),'%s','delimiter',',');
Param.Rov_detect = x{1}(:);clear x
%
str = 'rov param:';
a = '';
while ~contains(a,str)
    a = fgetl(fid);
    if a==-1;return;end
end
ii = find(a==',');
x = textscan(a(ii(1)+1:end),'%f','delimiter',',');
Param.Rov_param = x{1}(:);clear x
%
str = 'ref detect:';
a = '';
while ~contains(a,str)
    a = fgetl(fid);
    if a==-1;return;end
end
ii = find(a==',');
x = textscan(a(ii(1)+1:end),'%s','delimiter',',');
Param.Ref_detect = x{1}(:);clear x
%
str = 'ref param:';
a = '';
while ~contains(a,str)
    a = fgetl(fid);
    if a==-1;return;end
end
ii = find(a==',');
x = textscan(a(ii(1)+1:end),'%s','delimiter',',');
Param.Ref_param = x{1}(:);clear x
%
str = 'Low V ID:';
a = '';
while ~contains(a,str)
    a = fgetl(fid);
    if a==-1;return;end
end
ii = find(a==',');
x = textscan(a(ii(1)+1:end),'%f','delimiter',',');
Param.Low_voltage_ID = x{1}(:);clear x
%
str = 'Seg data len:';
a = '';
while ~contains(a,str)
    a = fgetl(fid);
    if a==-1;return;end
end
ii = find(a==',');
x = textscan(a(ii(1)+1:end),'%f','delimiter',',');
Param.Segment_data_length = x{1}(:);clear x
%
str = 'Exported seconds:';
a = '';
while ~contains(a,str)
    a = fgetl(fid);
    if a==-1;return;end
end
ii = find(a==',');
x = textscan(a(ii(1)+1:end),'%f','delimiter',',');
Param.Exported_seconds = x{1}(:);clear x
%
str = 'Sample rate:';
a = '';
while ~contains(a,str)
    a = fgetl(fid);
    if a==-1;return;end
end
ii = find(a==',');
x = textscan(a(ii(1)+1:end),'%f','delimiter',',');
Param.fs = x{1}(:);clear x
%
str = 'CFE P-P ';
a = '';
while ~contains(a,str)
    a = fgetl(fid);
    if a==-1;return;end
end
ii = find(a==',');
x = textscan(a(ii(1)+1:end),'%f','delimiter',',');
Param.CFE_sens = x{1}(:);clear x
%
str = 'CFE Width';
a = '';
while ~contains(a,str)
    a = fgetl(fid);
    if a==-1;return;end
end
ii = find(a==',');
x = textscan(a(ii(1)+1:end),'%f','delimiter',',');
Param.CFE_width = x{1}(:);clear x
%
str = 'CFE Refractory';
a = '';
while ~contains(a,str)
    a = fgetl(fid);
    if a==-1;return;end
end
ii = find(a==',');
x = textscan(a(ii(1)+1:end),'%f','delimiter',',');
Param.CFE_refractory = x{1}(:);clear x




%% Signals
fseek(fid,0,-1);
str = 'rov trace:';
a = '';
while ~contains(a,str)
    a = fgetl(fid);
    if a==-1;return;end
end
x = textscan(a,'%s','delimiter',',');
Labels = x{1}(2:end);

a = '';
while ~contains(a,str)
    a = fgetl(fid);
    if a==-1;return;end
end
celltraces = textscan(fid, repmat('%f',[1,Param.N+1]), 'delimiter', ',','collectoutput',1);
signals = celltraces{1}(:,2:end); % First column is ROV (empty)

% =
str = 'ref trace:';
a = '';
while ~contains(a,str)
    a = fgetl(fid);
    if a==-1;return;end
end
celltraces = textscan(fid, repmat('%f',[1,Param.N+1]), 'delimiter', ',','collectoutput',1);
Other_signals.signals_ref = celltraces{1}(:,2:end); % First column is ROV (empty);

str = 'spare1 trace:';
a = '';
while ~contains(a,str)
    a = fgetl(fid);
    if a==-1;return;end
end
celltraces = textscan(fid, repmat('%f',[1,Param.N+1]), 'delimiter', ',','collectoutput',1);
Other_signals.signals_spare1 = celltraces{1}(:,2:end); % First column is ROV (empty);

str = 'spare2 trace:';
a = '';
while ~contains(a,str)
    a = fgetl(fid);
    if a==-1;return;end
end
celltraces = textscan(fid, repmat('%f',[1,Param.N+1]), 'delimiter', ',','collectoutput',1);
Other_signals.signals_spare2 = celltraces{1}(:,2:end); % First column is ROV (empty);

% =
str = 'spare3 trace:';
a = '';
while ~contains(a,str)
    a = fgetl(fid);
    if a==-1;return;end
end
celltraces = textscan(fid, repmat('%f',[1,Param.N+1]), 'delimiter', ',','collectoutput',1);
Other_signals.signals_spare3 = celltraces{1}(:,2:end); % First column is ROV (empty);

% 
l = cell(1,Param.N);
for i = 1:length(l)
    l{i} = ['P',num2str(Param.Point_number(i)),'-',Labels{i}(7:end)];
end

Param.Label = l;
Param.xyz = xyz;
Param.xyz_surfP = xyz_surfP;
Param.Other_signals = Other_signals;

Param = orderfields(Param);
fclose all;

% % formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
% formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f';
% textscan(fid, '%[^\n\r]', 81-1, 'ReturnOnError', false);
% dataArray = textscan(fid, formatSpec,1076,'Delimiter', ',', 'EmptyValue' ,NaN,'ReturnOnError', false);
%
%
% signals =fscanf(fid,'%f');


%
% str_target1 = 'Wave Names and ( x y z ) coordinates ';
% str_target2 = 'Begin data';
%
% Header = fgetl(fid);
% Stduy = fgetl(fid);
% a = fgetl(fid);
% jj = find(a==':');
% Time_start = a(jj(1)+2:end);
% a = fgetl(fid);
% Time_end = a(jj(1)+2:end);
%
% jj1 = find(Time_start==':');
% x1 = str2double(Time_start(1:jj1(1)-1))*(60*60) + str2double(Time_start(jj1(1)+1:jj1(2)-1))*(60) + str2double(Time_start(jj1(2)+1:jj1(3)-1)) + str2double(Time_start(jj1(3)+1:end))/1000;
% jj1 = find(Time_end==':');
% x2 = str2double(Time_end(1:jj1(1)-1))*(60*60) + str2double(Time_end(jj1(1)+1:jj1(2)-1))*(60) + str2double(Time_end(jj1(2)+1:jj1(3)-1)) + str2double(Time_end(jj1(3)+1:end))/1000;
% Time_Duration = x2-x1; % seconds
% clear jj* x2 x1
%
% a = fgetl(fid);
% a = fgetl(fid);
% jj = find(a==':');
% Channels_num = str2double(deblank(a(jj+1:end)));
% a = fgetl(fid);
% jj = find(a==':');
% Samples_num = str2double(deblank(a(jj+1:end)));
%
%
% c=1;
% while c
%     a=fgetl(fid);
%     if ~isempty(strfind(a, str_target1))
%         c = 0;
%     end
%     if a==-1
%        error('Ensite Version may be unsupported')
%     end
% end
%
% a=fgetl(fid);
% ind = 0;indv=0;
% Label = cell(Channels_num,1);
% xyz = nan(1,3);
% while isempty(strfind(a, str_target2))
%     ind = ind+1;
%
%     a=fgetl(fid);
%     iv = findstr(a,'Virtual');
%
%     if isempty(iv)
%         jj = find(a=='=');
%         Label{ind} = a(jj+1:end);
%     else
%         indv = indv+1;
%
%         Label{ind} = ['Virtual-',num2str(indv,'%1.2d')];
%         jj = find(a=='(');
%         x = textscan(a(jj(2)+1:end-1),'%f%f%f','CollectOutput',1);
%         xyz(indv,:) = x{1};
%     end
% end
% Label(cellfun(@isempty,Label)) = [];
%
%
% a=fgetl(fid);
% signals =fscanf(fid,'%f',[Channels_num,inf]);
% signals = signals';
% % t = linspace(0,Time_Duration,Samples_num);
% % fs =  1/nanmean(diff(t));
% fs =  1200; % Hz
% t = [0:size(signals,1)-1]/fs;
% % reorder Labes
% [Label_sort,iis] = sort(Label);
% signals = signals(:,iis);
%
% Param.frequency = fs;
% Param.t = t;
% Param.Header = Header;
% Param.Stduy = Stduy;
% Param.Time_start = Time_start;
% Param.Time_end = Time_end;
% Param.Time_Duration = Time_Duration;
% Param.Label = Label_sort;
% Param.Virtual_xyz = xyz;
%
% fclose all;


function [signals,Param] = load_DxL_Precision(filename);

% Example:
% filename = 'G:\Others\Velocity_Export\study_dwsG600573_2018_02_13_13_10_26\2018_03_16_17_57_16\DxL_14.csv';
% [signals,Param] = load_DxL_Precision(filename);
Param = [];
fid=fopen(filename);


str = 'Total number of data points (columns)';
while 1
    a = fgetl(fid);
    if ~isempty(strfind(a,str));break;end
    if a==-1;return;end
end
ii = find(a==',');
x = deblank(a(ii(end)+1:end));
N = str2double(x);
%
str = 'pt number:';
while 1
    a = fgetl(fid);
    if ~isempty(strfind(a,str));break;end
    if a==-1;return;end
end
ii = find(a==',');
x = textscan(a(ii(1)+1:end),'%f','delimiter',',','collectoutput',1);
Point_number = x{1};clear x
%
xyz = nan(N,3);
str = 'roving x:';
while 1
    a = fgetl(fid);
    if ~isempty(strfind(a,str));break;end
    if a==-1;return;end
end
ii = find(a==',');
x = textscan(a(ii(1)+1:end),'%f','delimiter',',','collectoutput',1);
xyz(:,1) = x{1}(:);clear x
str = 'roving y:';
while 1
    a = fgetl(fid);
    if ~isempty(strfind(a,str));break;end
    if a==-1;return;end
end
ii = find(a==',');
x = textscan(a(ii(1)+1:end),'%f','delimiter',',');
xyz(:,2) = x{1}(:);clear x
str = 'roving z:';
while 1
    a = fgetl(fid);
    if ~isempty(strfind(a,str));break;end
    if a==-1;return;end
end
ii = find(a==',');
x = textscan(a(ii(1)+1:end),'%f','delimiter',',');
xyz(:,3) = x{1}(:);clear x
% =
xyz_surfP = nan(N,3);
str = 'surfPt x:';
while 1
    a = fgetl(fid);
    if ~isempty(strfind(a,str));break;end
    if a==-1;return;end
end
ii = find(a==',');
x = textscan(a(ii(1)+1:end),'%f','delimiter',',');
xyz_surfP(:,1) = x{1}(:);clear x
str = 'surfPt y:';
while 1
    a = fgetl(fid);
    if ~isempty(strfind(a,str));break;end
    if a==-1;return;end
end
ii = find(a==',');
x = textscan(a(ii(1)+1:end),'%f','delimiter',',');
xyz_surfP(:,2) = x{1}(:);clear x
str = 'surfPt z:';
while 1
    a = fgetl(fid);
    if ~isempty(strfind(a,str));break;end
    if a==-1;return;end
end
ii = find(a==',');
x = textscan(a(ii(1)+1:end),'%f','delimiter',',');
xyz_surfP(:,3) = x{1}(:);clear x
% - 
str = 'utilized:';
while 1
    a = fgetl(fid);
    if ~isempty(strfind(a,str));break;end
    if a==-1;return;end
end
ii = find(a==',');
x = textscan(a(ii(1)+1:end),'%f','delimiter',',');
Utilized = x{1}(:);clear x


% = 
% - Signal length (s)
str = 'Seg data len';
while 1
    a = fgetl(fid);
    if ~isempty(strfind(a,str));break;end
    if a==-1;return;end
end
Param.Segment_data_length = str2double(a(find(a==',')+1:end));
% - Export length (s)
str = 'Exported seconds:';
while 1
    a = fgetl(fid);
    if ~isempty(strfind(a,str));break;end
    if a==-1;return;end
end
Param.Exported_length = str2double(a(find(a==',')+1:end));
% - Frequency (s)
str = 'Sample rate:';
while 1
    a = fgetl(fid);
    if ~isempty(strfind(a,str));break;end
    if a==-1;return;end
end
Param.Frequency = str2double(a(find(a==',')+1:end));
% -
Param.Exported_length = str2double(a(find(a==',')+1:end));

%%
% - Frequency (s)
fseek(fid,0,-1);
str = 'ref LAT:';
while 1
    a = fgetl(fid);
    if ~isempty(strfind(a,str));break;end
    if a==-1;return;end
end
x = textscan(a,'%s%f','delimiter',',');
Param.Index_ref = x{2}(2:end);
a = fgetl(fid);
x = textscan(a,'%s%f','delimiter',',');
Param.Index_rov = x{2}(2:end);
%
fseek(fid,0,-1);
str = 'rov detect:';
while 1
    a = fgetl(fid);
    if ~isempty(strfind(a,str));break;end
    if a==-1;return;end
end
x = textscan(a,'%s','delimiter',',');
Param.Index_type = x{1}(2:end);






% = signals
fseek(fid,0,-1);
str = 'rov trace:';
while 1
    a = fgetl(fid);
    if ~isempty(strfind(a,str));break;end
    if a==-1;return;end
end
x = textscan(a,'%s','delimiter',',');
Labels = x{1}(2:end);
while 1
    a = fgetl(fid);
    if ~isempty(strfind(a,str));break;end
    if a==-1;return;end
end
celltraces = textscan(fid, repmat('%f',[1,N+1]), 'delimiter', ',','collectoutput',1);
signals = celltraces{1}(:,2:end); % First column is ROV (empty)

% =
str = 'ref trace:';
while 1
    a = fgetl(fid);
    if ~isempty(strfind(a,str));break;end
    if a==-1;return;end
end
celltraces = textscan(fid, repmat('%f',[1,N+1]), 'delimiter', ',','collectoutput',1);
Other_signals.signals_ref = celltraces{1}(:,2:end); % First column is ROV (empty);

str = 'spare1 trace:';
while 1
    a = fgetl(fid);
    if ~isempty(strfind(a,str));break;end
    if a==-1;return;end
end
celltraces = textscan(fid, repmat('%f',[1,N+1]), 'delimiter', ',','collectoutput',1);
Other_signals.signals_spare1 = celltraces{1}(:,2:end); % First column is ROV (empty);

str = 'spare2 trace:';
while 1
    a = fgetl(fid);
    if ~isempty(strfind(a,str));break;end
    if a==-1;return;end
end
celltraces = textscan(fid, repmat('%f',[1,N+1]), 'delimiter', ',','collectoutput',1);
Other_signals.signals_spare2 = celltraces{1}(:,2:end); % First column is ROV (empty);

% =
str = 'spare3 trace:';
while 1
    a = fgetl(fid);
    if ~isempty(strfind(a,str));break;end
    if a==-1;return;end
end
celltraces = textscan(fid, repmat('%f',[1,N+1]), 'delimiter', ',','collectoutput',1);
Other_signals.signals_spare3 = celltraces{1}(:,2:end); % First column is ROV (empty);

% 
l = cell(1,length(Point_number));
for i = 1:length(Point_number)
    l{i} = ['P',num2str(Point_number(i)),'-',Labels{i}(7:end)];
end

Param.Label = l;
Param.xyz = xyz;
Param.xyz_surfP = xyz_surfP;
Param.Other_signals = Other_signals;
Param.Point_number = Point_number;
Param.Utilized = Utilized;
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


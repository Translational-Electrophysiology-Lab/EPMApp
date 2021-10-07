%% Function to export data from ENSITE file.
% Only suitable for St. Jude Medical. File Revision : Velocity 2.0
%% MO 04/2014
function [Data,geo_virtuals]=load_data_from_ensite_mo_v2(filename1)

% [Data]=load_data_from_ensite_mo(filename1)
% INPUT
%   - filenemame1: string with path+name+format
% OUTPUT
%   - Data: strunct array with fields:
%       Type : type of data (ECG, EGMs, Virtual)
%       wavenames : name of waves
%       traces : data (as many columns as length(wavenames))
%       Times : matrices with time info
%       Time_name : type of time variable
%%

hb = waitbar(0,'Importing Data ...');
fid=fopen(filename1);

%%
geo_virtuals.Labels = [];
geo_virtuals.xyz = [];

%% open file. check if this is ENSITE file. get format type
% i.e. find "St. Jude Medical Virtual  data export; file format revision 3"
% in first row of the file


%%
seekstr = 't_dws,t_secs,t_usecs,t_ref';  %channel names header line - next line start of the actual data
a=fgetl(fid);
Data.Type = a(findstr(a,':')+2:end);
while 1
    a=fgetl(fid);
    if ~ischar(a), display('Broken file'); return; , end % since we dont know number of electrodes halt the function!
    if (strfind(a, seekstr)), break, end
end

%%
a(end+1) = ',';
wavenames_temp = line2cellArray(a); % Name of waveform
dataposition = ftell(fid);
tempdata = line2cellArray(fgetl(fid)); % read the first line of data to get the number of columns


%go back to the data position (end of the header
fseek(fid,dataposition,'bof');
%makeformat string, put %s as number of data channels
formatstr = '%*s '; % skip the first one
ind=0;
Data.wavenames = cell(1);
for kk=2:length(tempdata)
    if ~isempty(regexp(wavenames_temp{kk},'_ds'))|~isempty(regexp(wavenames_temp{kk},'_ps'))
        formatstr = [formatstr,'%*f '];
    else
        ind=ind+1;
        Data.wavenames{ind} = wavenames_temp{kk};
        formatstr = [formatstr,'%f '];
    end
end
waitbar(1/3,hb)
%read all data into the one large cellarray
celltraces = textscan(fid, formatstr, 'delimiter', ',','collectoutput',1);
%     wavenames = wavenames(1:size(celltraces,2));
%     ar(1/3,hb);
%     iiok = cellfun(@isnumeric,celltraces);
%     traces = cell2mat(celltraces(iiok));
waitbar(2/3,hb)
%     wavenames = wavenames(2:length(celltraces));
%     frequency = 2000;%round(1/(traces(2,1) - traces(1,1)));

%%
Data.Times = celltraces{1}(:,strncmp(Data.wavenames,'t_',2));%clear celltraces
Data.Times_name = Data.wavenames(strncmp(Data.wavenames,'t_',2));%clear celltraces
Data.traces = celltraces{1}(:,~strncmp(Data.wavenames,'t_',2));%clear celltraces
Data.wavenames(strncmp(Data.wavenames,'t_',2))=[];%clear celltraces
Data.traces(end,:)=[]; % remove nan due to 'eof'
Data.Times(end,:)=[]; % remove nan due to 'eof'
%     wavename2 = wavenames(~strncmp(wavenames,'t_',2));
%     traces2 = traces2(:,~ [~cellfun(@isempty,regexp(wavename2,'_ds')) | ~cellfun(@isempty,regexp(wavename2,'_ps'))]);
%     wavename2 = wavename2(~ [~cellfun(@isempty,regexp(wavename2,'_ds')) | ~cellfun(@isempty,regexp(wavename2,'_ps'))]);

%     Data.traces = traces2;
%     Data.wavename = wavename2;
waitbar(1,hb)

if isequal(Data.Type,'EPcathBIO_RAW')
    fseek(fid,0,'bof');
    seekstr1 = 'used channels';  % mo
    
    while 1
        a=fgetl(fid);  % "fgetLINE" not "fgetONE"
        if ~ischar(a), display('Broken file'); return; , end % since we dont know number of electrodes halt the function!
        if (strfind(a, seekstr1)), break, end
    end
    ind = 0;
    while ~isequal(a,'')&ind<100
        ind = ind+1;
        a = fgetl(fid);
        Data.Info_channels{ind} = a;
    end
end

%% xyz for grid
fseek(fid,0,-1)
a=fgetl(fid);
seekstr = 'Virtual locations';  %channel names header line - next line start of the actual data
a=fgetl(fid);
Data.Type = a(findstr(a,':')+2:end);
ind = 0;
while ind<200
    ind = ind+1;
    a=fgetl(fid);
    if (strfind(a, seekstr)), break, end
end
if ind<200
    points = textscan(fid,'%f%f%f', 'delimiter', ',','collectoutput',1);
    xyz = points{1};
    geo_virtuals.Labels = Data.wavenames;
    geo_virtuals.xyz = xyz;
else
    fseek(fid,0,-1)
end

%% xyz for virtuals
if isequal(Data.Type,'Virtuals')
    fseek(fid,0,-1)
    a=fgetl(fid);
    seekstr = 'User placed virtual';  %channel names header line - next line start of the actual data
    ind = 0;
    while ind<200
        ind = ind+1;
        a=fgetl(fid);
        if (strfind(a, seekstr)), break, end
    end
    if ind<200
        xx = textscan(fid,'%s%f%f%f', 'delimiter', ',','collectoutput',1);
        iiok = ~isnan(sum(xx{2},2));
        geo_virtuals.Labels = xx{1}(iiok);
        geo_virtuals.xyz = xx{2}(iiok,:);
    else
        fseek(fid,0,-1)
    end
end

%% Close the file
fclose(fid);

%now a function to help clean up any NaN's in the file - as with velocity
%with different (?) export function - for Justine;

disp('Check import')
close(hb)


%% some help functions
function output=line2cellArray(line) %create

nchannels=strfind(line,',');
previndex = 1;

for kk=1:length(nchannels)
    output{kk} = line(previndex:nchannels(kk)-1);
    previndex = nchannels(kk)+1;
end

function output=line2array(line)

nchannels=strfind(line,',');
previndex = 1;

for kk=1:length(nchannels)
    if kk == 1
        output(kk) = convert2seconds(line(previndex:nchannels(kk)-1));
    else
        output(kk) = str2double(line(previndex:nchannels(kk)-1));
    end
    previndex = nchannels(kk)+1;
end

function finaloutput=convert2seconds(tvector)

% Example
% With "jump"                   Without jump
% 22.42.05.655,4.67923
% 22.42.05.656,4.67972  |    22.42.05.656,4.67972
% 22.42.05.661,0        |    22.42.05.656,4.68021
% 22.42.05.662,0.000491 |    22.42.05.657,4.6807

finaloutput = [];
for i=1:length(tvector)
    text1 = tvector{i};
    ndots = strfind(text1,'.');
    output = str2double(text1(1:ndots(1)-1))*60*60;
    output = output + str2double(text1(ndots(1)+1:ndots(2)-1))*60;
    output = output + str2double(text1(ndots(2)+1:end));
    finaloutput = [finaloutput;output];
end


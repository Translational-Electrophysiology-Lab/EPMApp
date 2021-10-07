function [signals,Param] = LoadBardSignals_moCathLab(filename,do_only_header);
% function [signals,valid_exp,frequency,start_time,end_time,scale] = LoadBardSignals_moCathLab(filename,bard_num,geo_fname,chan_fname);
% The same as LoadBardSignals.m except for calling outputs instead of
% saving them
%%
%save [signals, exp_channel, exp_channel_elect_name, exp_channel_index, frequency, start_time, end_time] into sig_ori_fname
%input file has bard_num channels
% bard_num = 1: Bard 1, export 160 channels
% bard_num = 2: Bard 1, export 80 channels
%global sig_ori_fname
%%

% outputfiles = 'output/outputfilenames.mat';
% load(outputfiles,'sig_ori_fname','chan_fname','geo_fname');
% sig_ori_fname
% chan_fname
% geo_fname

if nargin<2
    do_only_header = 0;
end

nchannels = [];  % total num of channels
Param.frequency = [];
Param.scale = [];
Param.signals = [];
Param.range = [];
Param.Label = cell(0,1);
Param.units = cell(0,1);
Param.Filter_BW = cell(0,2);
Param.filename_ori = {filename};
fid=fopen(filename);

seekstr = '[Data]';
while 1
    a=fgetl(fid);
    if ~ischar(a), display('Broken file'); return; end
    if (strfind(a, seekstr)), break, end
end
dataposition = ftell(fid);

fseek(fid,0,'bof');
while ftell(fid) < dataposition
    a = fgetl(fid);
    if (strfind(a,'Channels exported'))
        tmp = strfind(a,':');
        nchannels = round(str2num(a(tmp+2:end)));
        Param.Chann_Exp_Number = nchannels;
        fprintf('Exported %i channels ',nchannels)
    end
    if(strfind(a,'Start time: '))
        tmp = strfind(a,':');
        fprintf('from %s', a(tmp+2:end));
        tmp = sscanf(a(tmp+2:end),'%f:%f:%f');
        Param.start_time_string = a;
        Param.start_time = tmp(1)*3600+tmp(2)*60+tmp(3); %in sec
    end
    if(strfind(a,'End time: '))
        tmp = strfind(a,':');
        fprintf(' to %s.\n', a(tmp+2:end));
        tmp = sscanf(a(tmp+2:end),'%f:%f:%f');
        Param.end_time = tmp(1)*3600+tmp(2)*60+tmp(3); % in sec
        Param.end_time_string = a; % in sec
    end
    
    if strfind(a,'Data Format: ')
        tmp = strfind(a,':');
        Param.DataFormat = str2num(a(tmp+2:end-2));
    end
    
    if strfind(a,'Sample Rate: ')
        tmp = strfind(a,':');
        Param.frequency = str2num(a(tmp+2:end-2));
    end
    if strfind(a,'Range')
        tmp = strfind(a,':');
        y = a(tmp+2 : end);
        iy = 0;
        if ~isempty(deblank(y))
            while isempty(str2num(y(1:length(y)-iy)));
                iy = iy+1;
            end
            Param.range(end+1) = str2num(y(1:length(y)-iy));
            Param.units{end+1} = deblank(y(end-iy+1:end));
        else
            Param.range(end+1) = 5000;
            Param.units{end+1} = 'mV';
        end
    end
    if strfind(a,'Scale')
        tmp = strfind(a,':');
        Param.scale(end+1) = str2num(a(tmp+2:end));
    end
    if strfind(a,'Label')
        tmp = strfind(a,':');
        Param.Label(end+1) = {a(tmp+2:end)};
    end
     if strfind(a,'Low')
        tmp = strfind(a,':');
        Param.Filter_BW(end+1,1) = {a(tmp+2:end)};
        a = fgetl(fid);
        tmp = strfind(a,':');
        Param.Filter_BW(end,2) = {a(tmp+2:end)};
     end
     
end
fprintf('Total exported duration: %6.3f (sec).\n',Param.end_time-Param.start_time)

if do_only_header
    signals = [];
    return
end
fseek(fid,dataposition,'bof');
formatstr = '';
for i=1:nchannels
    formatstr = [formatstr,'%f '];
end

celltraces = textscan(fid,formatstr,'delimiter',',');
L = cellfun(@length,celltraces);
if sum(abs(diff(L)))==0
    signals = cell2mat(celltraces);% xiao
else
    Lm = min(L);
    A = cellfun(@(x) x(1:Lm),celltraces,'UniformOutput', false);
    signals = cell2mat(A);% xiao
    warning(['Channels have different length'])
    warning([' - Max different length = ',num2str(range(L))])
    
end


for i=1:nchannels
    %     signals(:,i) = power(10,Param.scale(i))*signals(:,i)*1000; %mV
    %     signals(:,i) = signals(:,i)*2^(Param.scale(i)); %mV 18/11/2016
    signals(:,i) = signals(:,i)/(2^15)*(Param.range(i)); %mV 18/11/2016
end
fclose(fid);



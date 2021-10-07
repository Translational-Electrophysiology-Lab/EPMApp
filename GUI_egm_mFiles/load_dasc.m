function [signals,ParamSig] = load_dasc(filename); 
% [signals,ParamSig] = load_dasc(filename,N);
% Read DASC files
% - Input: filename: String including path and format
%        : N: Scalar, indicating the number of channels recorded
% - Output: signals: Signals
%         : ParamSig.date = Date of recording
%         : ParamSig.time = Time of recording
%         : ParamSig.frequency = Sampling freq
%         : ParamSig.Label = Labels
%% rewrite the file and replace comma with points

file    = memmapfile( filename, 'writable', true );
comma   = uint8(',');
point   = uint8('.');
file.Data( transpose( file.Data==comma) ) = point;

% Open file
fid = fopen(filename);

strseek = 'Date and';
move_on = 1;
while move_on
    a = fgetl(fid);
    if ~isempty(strfind(a,strseek))
        move_on = 0;
    end
end
ii = find(a=='-');
ParamSig.date = a(ii(1)-2:ii(2)+4);
ii = find(a==':');
ParamSig.time = a(ii(1)-2:ii(2)+2);

%strseek = 'Type';
strseek = 'Time (s)';
move_on = 1;
while move_on
    a = fgetl(fid);
    if ~isempty(strfind(a,strseek))
        move_on = 0;
    end
end
a = fgetl(fid);
N = sum(isspace(a))-1;
clear a
clear ii

h = textscan(fid,['%s',repmat('%f',[1 N])],'collectoutput',1);

signals = h{2};

timev = h{1};
ii1 = strfind(timev,':');
ii2 = strfind(timev,'.');
t = nan(length(timev),1);
for i = 1:length(timev)
    th = 3600*str2double(timev{i}(ii1{i}(1)-2 : ii1{i}(1)-1));
    tm = 60*str2double(timev{i}(ii1{i}(2)-2 : ii1{i}(2)-1));
    ts = str2double(timev{i}(ii1{i}(2)+1 : ii1{i}(2)+2));
    tms = str2double(timev{i}(ii2{i}(1)+1 : ii2{i}(1)+3));
    
    t(i) = th+tm+ts+tms;
end

ParamSig.t = t;
ParamSig.frequency = 1000/round(nanmedian(diff(t)));
% ParamSig.Label = num2cell(1:N);



function [signals,Labels,Info] = LoadWaveformsFromCarto(loadname);

% dir_load = 'E:\UCL\Data&People\Data_Barts_VT\151202\Patient 2015_12_02(1)\Study 1\Export_Study-1-12_02_2015-21-07-41\';
% filename = '1-1-ReMap_P49_ECG_Export.txt';
% loadname = [dir_load,filename];

delimiter = ' ';
startRow = 1;
fid=fopen(loadname);
Info.ExportType = fgetl(fid);
Info.Gain = fgetl(fid);
Info.Title = fgetl(fid);

% Labels = textscan(fid, '%s', 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'EmptyValue' ,NaN, 'ReturnOnError', false);
% Labels = textscan(fid,'%s',50,'MultipleDelimsAsOne', true,'Delimiter', ' ')%;,'endofline','\n');
L = fgetl(fid);
a = find(L~=' ');
ii = [0 find(diff(a)>1.1) length(a)];
Labels = cell(1,length(ii)-1);
for i = 1:length(ii)-1
    Labels(i) = {L(a(ii(i)+1):a(ii(i+1)))};
end
formatSpec = [repmat('%f',[1,length(Labels)]),'%[^\n\r]'];
% formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
% Labels = textscan(fid,'%s','delimiter',',');
% Labels = fscanf(fid,'%s',[5 1]);
dataArray = textscan(fid, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false);
signals = [dataArray{1:end-1}];

%
ii = find(Info.Gain=='=');
gain = str2double(Info.Gain(ii+1:end));
Info.Gain_value = gain;
if ~isnan(gain)
    signals = signals*gain;
end
fclose(fid);

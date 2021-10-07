function [signals,ParamSig] = LoadAbbottsCSV(loadname);

fid=fopen(loadname);
a = fgetl(fid);
a = fgetl(fid);
ParamSig.frequency = str2double(a(14:end));
a = fgetl(fid);
ParamSig.resolution = str2double(a(28:end));
a = fgetl(fid);
N = str2double(a(33:end));

seekstr = 'Sample Number';
a=fgetl(fid);
while 1
    a=fgetl(fid);
    if ~ischar(a), display('Broken file'); return; , end % since we dont know number of electrodes halt the function!
    if (strfind(a, seekstr)), break, end
end
% a=fgetl(fid);

ii = find(a==',');
lab = cell(1,length(ii)+1);
lab{1} = a(1:ii(1)-1);
for i = 1:length(ii)-1
    lab{i+1} = a(ii(i)+1:ii(i+1)-1);
end
lab{end} = a(ii(end)+1:end);

% X = textscan(fid,repmat('%f',[1 N]),'CollectOutput',1);
N = sum(a==',')+1;
formatSpec = repmat('%s',[1 N]);

X = textscan(fid, formatSpec, 'Delimiter', ',', 'MultipleDelimsAsOne', true, 'EmptyValue' ,NaN,'CollectOutput',1);

signals = nan(size(X{1}));
hb = waitbar(0,'... Importing data');
for i = 1:size(signals,2)
    try
        signals(:,i) = str2double(X{1}(:,i));
    catch me
    end
    waitbar(i/size(signals,2));
end
close(hb);

ParamSig.t = signals(:,1)/ParamSig.frequency;
ParamSig.Label = lab(2:end);
signals(:,1) = [];
clear X

fclose(fid);



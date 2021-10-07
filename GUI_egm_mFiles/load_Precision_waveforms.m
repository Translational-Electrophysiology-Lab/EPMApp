function [signals,Param] = load_Precision_waveforms(filename);

fid=fopen(filename);

str = 'Export Data Element ';
while 1
    a = fgetl(fid);
    if contains(a,str);break;end
    if a==-1;return;end
end
ii = find(a==':');
x = deblank(a(ii(end)+2:end));
Param.Export_Data_Element = x;

fseek(fid,0,-1);
str = 'Export Data Element';
while 1
    a = fgetl(fid);
    if contains(a,str);break;end
    if a==-1;return;end
end
ii = find(a==':');
x = deblank(a(ii(end)+2:end));
Param.Export_Data_Element = x;

fseek(fid,0,-1);
str = 'Export from Study ';
while 1
    a = fgetl(fid);
    if contains(a,str);break;end
    if a==-1;return;end
end
ii = find(a==':');
x = deblank(a(ii(end)+2:end));
Param.Export_from_Study = x;

fseek(fid,0,-1);
str = 'Export files Stored in Dir';
while 1
    a = fgetl(fid);
    if contains(a,str);break;end
    if a==-1;return;end
end
ii = find(a==':');
x = deblank(a(ii(end)+2:end));
Param.Dir_ori = x;
% -
fseek(fid,0,-1);
str = 'Export Start Time (h:m:s.msec) : ';
while 1
    a = fgetl(fid);
    if contains(a,str);break;end
    if a==-1;return;end
end
ii = find(a==')');
x = deblank(a(ii(end)+4:end));
Param.Time_Start = x;

% Highpass filter
while 1
    if contains(a,'Highpass');
        ii = find(a==':');
        Param.Highpass = a(ii+2:end);
        break;
    end
    if contains(a,'Number of samples ');
        Param.Highpass = 'None';
        break
    end
    a = fgetl(fid);
end
% Lowpass filter
while 1
    if contains(a,'Lowpass');
        ii = find(a==':');
        Param.Lowpass = a(ii+2:end);
        break;
    end
    if contains(a,'Number of samples ');
        Param.Lowpass = 'None';
        break
    end
    a = fgetl(fid);
    
end
% Notch filter
while 1
    if contains(a,'Notch');
        ii = find(a==':');
        Param.Notch = a(ii+2:end);
        break;
    end
    if contains(a,'Number of samples ');
        Param.Notch = 'None';
        break
    end
    a = fgetl(fid);
end

% Number of catheters
while 1
    if contains(a,'Number of Catheters');
        ii = find(a==':');
        Param.N_catheters = a(ii+2:end);
        break;
    end
    if contains(a,'Number of samples ')
        Param.N_catheters = 'None';
        break
    end
    a = fgetl(fid);
end

% -
fseek(fid,0,'bof');
str = 'used channels : ';
do_ECG=0;
while 1
    a = fgetl(fid);
    if contains(a,str);break;end
    if contains(a,'Number of samples ');
        do_ECG=1;
        ii = find(a==',');
        Param.Samples = str2double(a(ii+1:end));
        break;
    end
end

if ~do_ECG
    a = fgetl(fid);
    a = fgetl(fid);
    A =[];
    while ~isempty(a)
        x = textscan(a,'%s','delimiter',',','collectoutput',1);
        A = [A;x{1}'];
        a = fgetl(fid);
    end
    Param.Label = cell(1,size(A,1));
    chan = cell(1,size(A,1));
    for i = 1:length(A)
        Param.Label{i} = [A{i,2},'-',A{i,3}];
        chan{i} = ['c',A{i,1}];
    end
    
    % used bipol channels :
    fseek(fid,0,'bof');
    str = 'used bipol channels : ';
    while 1
        a = fgetl(fid);
        if contains(a,str)
            do_bipol = 1;
            break;
        end
        if contains(a,'Number of samples ')
            do_bipol = 0;
            break;
        end
    end
    if do_bipol
        a = fgetl(fid);
        B = textscan(fid,'%n %n %n', 'delimiter', ',','collectoutput',1);
        chan = cell(1,size(B{1},1));
        for i = 1:size(B{1},1)
            chan{i} = ['c',num2str(B{1}(i,1))];
        end
    end
    
    fseek(fid,0,'bof');
    str = 'Number of samples ';
    while 1
        a = fgetl(fid);
        if contains(a,str);break;end
        if a==-1;return;end
    end
    ii = find(a==',');
    Param.Samples = str2double(a(ii+1:end));
    a = fgetl(fid);
    
    % - Signals
    x = textscan(a,'%s','delimiter',',','collectoutput',0);
    if do_bipol
        ii2 = nan(1,length(chan));
        Param.Label = cell(1,length(chan));
        for i = 1:length(ii2)
            ii2(i) = find(contains(x{1},[chan{i},'_ds']))-1;
            d = x{1}{ii2(i)};ii =find(d=='_');
            Param.Label{i} = d(1:ii-1);
        end
    else
        [~,~,ii2] = intersect(chan,x{1},'stable');
    end
    
    R = repmat({'%*s '},[1,length(x{1})]);
    R(ii2) = {'%n'};
else
    a = fgetl(fid);
    x = textscan(a,'%s','delimiter',',','collectoutput',0);
    ii2 = find(contains(x{1},'_ds'))-1;
    Param.Label = x{1}(ii2);
    R = repmat({'%*s '},[1,length(x{1})]);
    R(ii2) = {'%n'};
end
celltraces = textscan(fid,[R{:}], 'delimiter', ',','collectoutput',1);
signals = celltraces{1};
if size(signals,1)>Param.Samples
    signals(Param.Samples+1,:) = [];
end
fclose(fid);




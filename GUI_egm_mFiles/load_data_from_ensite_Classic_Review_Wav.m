function [signals,Param] = load_data_from_ensite_Classic_Review_Wav(filename1);

fid=fopen(filename1);
str_target1 = 'Wave Names and ( x y z ) coordinates ';
str_target2 = 'Begin data';

Header = fgetl(fid);
Stduy = fgetl(fid);
a = fgetl(fid);
jj = find(a==':');
Time_start = a(jj(1)+2:end);
a = fgetl(fid);
Time_end = a(jj(1)+2:end);

jj1 = find(Time_start==':');
x1 = str2double(Time_start(1:jj1(1)-1))*(60*60) + str2double(Time_start(jj1(1)+1:jj1(2)-1))*(60) + str2double(Time_start(jj1(2)+1:jj1(3)-1)) + str2double(Time_start(jj1(3)+1:end))/1000;
jj1 = find(Time_end==':');
x2 = str2double(Time_end(1:jj1(1)-1))*(60*60) + str2double(Time_end(jj1(1)+1:jj1(2)-1))*(60) + str2double(Time_end(jj1(2)+1:jj1(3)-1)) + str2double(Time_end(jj1(3)+1:end))/1000;
Time_Duration = x2-x1; % seconds
clear jj* x2 x1

a = fgetl(fid);
a = fgetl(fid);
jj = find(a==':');
Channels_num = str2double(deblank(a(jj+1:end)));
a = fgetl(fid);
jj = find(a==':');
Samples_num = str2double(deblank(a(jj+1:end)));


c=1;
while c
    a=fgetl(fid);
    if ~isempty(strfind(a, str_target1))
        c = 0;
    end
    if a==-1
       error('Ensite Version may be unsupported')
    end
end

a=fgetl(fid);
ind = 0;indv=0;
Label = cell(Channels_num,1);
xyz = nan(1,3);
while isempty(strfind(a, str_target2))
    ind = ind+1;
    
    a=fgetl(fid);
    iv = findstr(a,'Virtual');
    
    if isempty(iv)
        jj = find(a=='=');
        Label{ind} = a(jj+1:end);
    else
        indv = indv+1;
        
        Label{ind} = ['Virtual-',num2str(indv,'%1.2d')];
        jj = find(a=='(');
        x = textscan(a(jj(2)+1:end-1),'%f%f%f','CollectOutput',1);
        xyz(indv,:) = x{1};
    end
end
Label(cellfun(@isempty,Label)) = [];


a=fgetl(fid);
signals =fscanf(fid,'%f',[Channels_num,inf]);
signals = signals';
% t = linspace(0,Time_Duration,Samples_num);
% fs =  1/nanmean(diff(t));
fs =  1200; % Hz
t = [0:size(signals,1)-1]/fs;
% reorder Labes
[Label_sort,iis] = sort(Label);
signals = signals(:,iis);

Param.frequency = fs;
Param.t = t;
Param.Header = Header;
Param.Stduy = Stduy;
Param.Time_start = Time_start;
Param.Time_end = Time_end;
Param.Time_Duration = Time_Duration;
Param.Label = Label_sort;
Param.Virtual_xyz = xyz;

fclose all;


function Parameters = ReadProperties(fn)
A = h5read(fn,'/RawData/AcquisitionTaskDescription');

str_tot = {'SensitivityLowValue',...
    'SensitivityHighValue',...
    'SampleRate',...
    'Offset',...
    'NotchFilter',...
    'HighpassFilter',...
    'LowpassFilter',...
    'IsBipolar',...
    'IsTriggerChannel',...
    'DeviceNumber',...
    'DeviceName',...
    'LogicalChannelNumber',...
    'PhysicalChannelNumber'};

clear Parameters*
for j = 1:length(str_tot)
    str1 = ['<',str_tot{j},'>'];
    str2 = ['</',str_tot{j},'>'];
    ii1 = strfind(A{1},str1)+length(str1);
    ii2 = strfind(A{1},str2)-1;
    for k = 1:length(ii1)
        Parameters0.(str_tot{j})(k) = str2double(A{1}(ii1(k):ii2(k)));
    end
    if length(unique(Parameters0.(str_tot{j})))==1 | mean(isnan(Parameters0.(str_tot{j})))==1
    Parameters.(str_tot{j}) = Parameters0.(str_tot{j})(1);
    end
end


% str1 = '<SampleRate>';
% str2 = '</SampleRate>';
% ii1 = strfind(A{1},str1)+length(str1);
% ii2 = strfind(A{1},str2)-1;
% 
% for j = 1:length(ii1)
%     fs_tot(j) = str2double(A{1}(ii1(1):ii2(1)));
% end

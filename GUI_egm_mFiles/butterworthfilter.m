function x = butterworthfilter(ecg, Fs, Wn, high_or_low)

% x = butterworthfilter(ecg, Fs, Wn, high_or_low)
% ecg: vector
% Fs: frequency of the ecg signal
% Wn: if Wn = [high_pass_freq low_pass_freq], bandpass filter; otherwise,
% either high or low filter depending on high_or_low
% high_or_low: 'high', or, 'low'; string indicating whether it is high or
% low pass filter

if nargin < 3 || nargin > 4
    error('Wrong input arguments for butterworth filter\n',help(mfilename));
end

if nargin == 4
    if strcmp(high_or_low,'high')
        [b, a] = butter(3, Wn/(Fs/2), 'high');
    else
        [b, a] = butter(3, Wn/(Fs/2), 'low');
    end
    x = filtfilt(b,a,ecg);
end

if nargin == 3
    
    if Wn(1)==0
        [b, a] = butter(3,min(Wn(2)/(Fs/2),0.99), 'low');
        x = filtfilt(b,a,ecg);
        
    else
        if diff(Wn/(Fs/2))>0.01 % it seems that when bandpass is too small in comparison with Fs the filter gives numerical problems
            [b, a] = butter(3, Wn/(Fs/2),'bandpass');
            x = filtfilt(b,a,ecg);
        else
            [b, a] = butter(3, Wn(2)/(Fs/2),'low');
            x = filtfilt(b,a,ecg);
            [b, a] = butter(3, Wn(1)/(Fs/2),'high');
            x = filtfilt(b,a,x);
        end
    end
    
end




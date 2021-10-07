function [spikes] = find_pacing_spikes_mo(signals,Fs,sensthr,do_control,minCL)

if nargin<3
    sensthr = 0.75;
    minCL = 150;
end
if nargin<4
    do_control = 0;
    minCL = 150;    
end
if nargin<5
    minCL = 150; % delete spikes whihc correspond to a CL lower than minCL [ms]
end
% [spikes] = find_pacing_spikes_mo(signals,Fs)
% Michele Orini 2014
%%

spiky = signals;

spikydiff1 = diff(spiky);
spiky = diff(signals,3);
spiky = abs([spiky(1)*ones(3,1);spiky]);
spiky = filtfilt([1 1 1],3,spiky);
thresh = (nanmedian(spiky)*5 + sensthr*max(spiky))/2;
ll = find(sign(diff(spiky>thresh))==1);

ll(diff(ll)<minCL/1000*Fs)=[];
spikes1 = nan(1,length(ll));
spikes2 = nan(1,length(ll));
spikes3 = nan(1,length(ll));
spikes4 = nan(1,length(ll));
for i = 1:length(ll)
    H = round(ll(i)-10/1000*Fs) : round(ll(i)+10/1000);
    H(H<1)=[];
    [~,ii] = max(abs(spikydiff1(H)));
    spikes1(i) = ii+H(1)-1;
    
    [~,ii] = max(abs(signals(H)));
    spikes2(i) = ii+H(1)-1;
    
    [~,ii] = max(abs(spiky(H)));
    spikes3(i) = ii+H(1)-1;
    
    [~,ii] = max(spikydiff1(H));
    spikes4(i) = ii+H(1)-1;
end

d2 = diff(spikes2);
d2(d2 > nanmedian(d2)+20 | d2 < nanmedian(d2)-20)=[];
d1 = diff(spikes1);
d1(d1 > nanmedian(d1)+20 | d1 < nanmedian(d1)-20)=[];
d3 = diff(spikes3);
d3(d3 > nanmedian(d3)+20 | d3 < nanmedian(d3)-20)=[];
d4 = diff(spikes4);
d4(d4 > nanmedian(d4)+20 | d4 < nanmedian(d4)-20)=[];

[~,ii] = min([std(d1),std(d2),std(d3),std(d4)]);
if ii==1
    spikes = spikes1;
elseif ii==2
    spikes = spikes2;
elseif ii==3
    spikes = spikes3;
elseif ii==4
    spikes = spikes4;
end

if do_control
    x2 = diff(spikes) - median(diff(spikes));
    L = 9;
    ff = filter([ones(1,L)],L,x2);
    i0 = find(ff==0);i0(i0<L)=[];
    x2b = x2;
    ind = 1;
    while ind < length(i0)
        %     figure(1),
        %     plot(x2([i0(i)-L+1:i0(i)]));
        %     title([num2str(i0(i)),' ','sum=',num2str(sum(x2([i0(i)-L+1:i0(i)])))])
        %     pause
        x2b([i0(ind)-L+1:i0(ind)]) = 0 ;
        i0(i0 < i0(ind)+L) = [];
        ind = ind+1;
    end
    % ds2 = cumsum([spikes(2)-spikes(1) x2b]);
    ds2 = x2b + median(diff(spikes));
    spikes = cumsum([spikes(1) ds2]);
end

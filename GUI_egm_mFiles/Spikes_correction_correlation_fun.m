function [spikes_new,SNR_avg,Nabnormal_beats] = Spikes_correction_correlation_fun(sig,spikes,fs,Niter,CLmax,do_plot);
% [spikes_new,SNR_avg] = Spikes_correction_correlation_fun(sig,spikes,fs,Niter,CLmax);

do_cf3_normalization = 1;

if nargin<6
    do_plot = 0;
end
if nargin<5
    do_plot = 0;
    CLmax = 0.350; % ms
end
if nargin<4
    do_plot = 0;
    
    CLmax = 0.350; % ms
    Niter = 5;
end
%%

Amp_tot = linspace(0.1,0.90,Niter).^2;
% CLmax = 0.350;
%% remove baseline
sigf = medfilt1(sig(1:5:end,:),round(1.5*fs/5));
if rem(size(sig,1)-1,5)~=0
    sigf2 = interp1([1:5:size(sig,1) size(sig,1)]',[sigf;sigf(end,:)],[1:size(sig,1)]');
else
    sigf2 = interp1([1:5:size(sig,1)],[sigf],1:size(sig,1));
end
sig = sig-sigf2;
%%
spikes_samp = round(spikes/1000*fs);
maxQRS = round(CLmax*fs);
DeltaT = round(120/1000*fs);
sp= spikes_samp -DeltaT; %samples
Hend = DeltaT;
L=min([max(diff(sp)),maxQRS]) + Hend;
X = nan(L,length(sp),size(sig,2));
for i=1:length(sp)-1
    H = sp(i):sp(i+1)-1+ Hend;
    H(L+1:end)=[];H(H<1)=[];H(H>length(sig))=[];
    X(1:length(H),i,:) = sig(H,:);
end
H = sp(i+1):size(sig,1);
H(L:end)=[];
X(1:length(H),i+1,:) = sig(H,:);
clear sp

ii = [diff(spikes)>400 & diff(spikes)<1000];
Sm0 = squeeze(nanmean(X(:,ii,:),2));
Sm = nan(size(Sm0));
for j = 1:size(X,3)
    cc = corrcoef([Sm0(:,j),X(:,:,j)]);
    Sm(:,j) = nanmean(X(:,ii&cc(1,3:end)>0.70,j),2);
end

cf3_tot = nan(size(sig));
for j = 1:size(Sm,2)
    [c,tau] = xcorr(sig(:,j),Sm(:,j));
    c(1:find(tau==0)-DeltaT-1)=[];
    c(end-DeltaT+1:end) = [];
    c = c/nanstd(c);
    [B,A] = butter(3,50/fs);
    cf3 = filtfilt(B,A,c.^3);
    
    cf3 = cf3/(3*nanstd(cf3));
    
    cf3_tot(:,j) = cf3;
end

cf3 = nanmean(cf3_tot,2);
sp= spikes_samp -DeltaT; %samples

if do_cf3_normalization
    X2 = nan(size(X,1),size(X,2));
    for i=1:length(sp)-1
        H = sp(i):sp(i+1)-1+ Hend;
        H(L+1:end)=[];H(H<1)=[];H(H>length(sig))=[];
        X2(1:length(H),i) = cf3(H,:);
    end
    H = sp(i+1):size(sig,1);
    H(L:end)=[];
    X2(1:length(H),i+1) = cf3(H,:);
    clear sp
    
    mm = max(X2);
    mm = [abs(mm(1));mm(:);abs(mm(end))];
    mm = medfilt1(mm,10);
    xm = [1 spikes_samp size(sig,1)];
    mm2 = interp1(xm,mm,[1:size(sig,1)],'linear');
    cf3 = cf3(:)./abs(mm2(:));
    clear mm mm2 xm
end
cf3d = diff(cf3);cf3d2 =diff(cf3d);

spikes_samp_new_all= cell(1,length(Amp_tot));

for iAmp = 1:Niter
    clear iix*
    iix = find(cf3d(1:end-1).*cf3d(2:end)<0 & cf3d2<0 & cf3(3:end)> min([Amp_tot(iAmp),0.90]))+1;
    iix_tot_sort = unique(iix);
    ii = find(diff(iix_tot_sort)<10);
    iix_tot_sort2 = iix_tot_sort;
    for k = 1:length(ii)
        kk = find(abs(iix_tot_sort-iix_tot_sort(ii(k)))<50);
        iix_tot_sort2(kk) = round(nanmedian(iix_tot_sort(kk)));
    end
    iix_tot_sort2 = unique(iix_tot_sort2);
    rr_samp = diff(iix_tot_sort2);
    
    %         figure,
    %         ax(1)= subplot(211);
    %         plot(cf3),hold on,plot(iix_tot_sort2,cf3(iix_tot_sort2),'or')
    %         ax(2)= subplot(212);
    %         plot(iix_tot_sort2(2:end),diff(iix_tot_sort2),'-ok')
    %         linkaxes(ax,'x')
    
    RRmax = 1500;
    ii = find(abs(diff(rr_samp/fs*1000))>200 |rr_samp(2:end)/fs*1000<200|rr_samp(2:end)/fs*1000>2000); % ms
    iix_tot_sort3 = iix_tot_sort2;
    for k = 1:length(ii)
        xx = iix_tot_sort2(ii(k)) : iix_tot_sort2(ii(k)+1);
        cxx = cf3(xx);
        cxxd = diff(cxx);cxxd2 =diff(cxxd);
        
        hh = ii(k) + [-10 : 10];
        hh(hh<1|hh>length(iix_tot_sort2))=[];
        Alim = mean(cf3(iix_tot_sort2(hh)));
        rrdumm = diff(iix_tot_sort2(hh));
        rrdumm(rrdumm>RRmax*fs/1000 | rrdumm<300*fs/1000) = [];
        rrdumm = nanmedian(rrdumm);
        
        iix = find(cxxd(1:end-1).*cxxd(2:end)<0 & cxxd2<0 & cxx(3:end)> Alim*0.20)+1;
        %         iix = find(cxxd(1:end-1).*cxxd(2:end)<0 & cxxd2<0 & cxx(3:end)> 0.01)+1;
        
        iiko = diff([xx(1)  xx(iix)])>rrdumm*2 | diff([xx(1)  xx(iix)])<rrdumm*0.50;
        iix(iiko) = [];
        iix_tot_sort3 = sort([iix_tot_sort3;xx(iix(:))']);
        %         figure,plot(cf3_tot(:,j)),hold on,plot(xx,cf3_tot(xx,j),'r'),plot(xx(iix),cf3_tot(xx(iix),j),'ok')
        
    end
    
    
    %         figure,
    %         ax(1)= subplot(211);
    %         plot(cf3),hold on,plot(iix_tot_sort3,cf3(iix_tot_sort3),'or')
    %         ax(2)= subplot(212);
    %         plot(iix_tot_sort3(2:end),diff(iix_tot_sort3),'-ok')
    %         linkaxes(ax,'x')
    
    
    
    iix_tot_sort4 = iix_tot_sort3;
    T = round(350/1000*fs);
    ii = find(diff(iix_tot_sort4)<T);
    cf3m = cf3;
    % figure,plot(cf3),hold on,plot(iix_tot_sort3,cf3(iix_tot_sort3),'or')
    for k = 1:length(ii)
        kk = find(abs(iix_tot_sort4-iix_tot_sort4(ii(k)))<=T);
        [s,iks]= sort(cf3m(iix_tot_sort4(kk))/max(cf3m(iix_tot_sort4(kk))),'descend');
        if s(2)<0.66;
            iix_tot_sort4(kk) = iix_tot_sort4(kk(iks(1)));
        else
            iix_tot_sort4(kk) = round(nanmedian(iix_tot_sort4(kk)));
        end
    end
    iix_tot_sort4 = unique(iix_tot_sort4);
    rr_samp4 = diff(iix_tot_sort4);
    
    spikes_samp_new_all{iAmp} = iix_tot_sort4;
    % figure,plot(cf3),hold on,plot(iix_tot_sort4,cf3(iix_tot_sort4),'or')
    
end

clear iix_tot*
[m,iSD] = min(cellfun(@(x) nanstd(diff(x,2)),spikes_samp_new_all));
Nabnormal_beats = cellfun(@(x) sum(abs(diff(x,2))>250),spikes_samp_new_all);
Nabnormal_beats = Nabnormal_beats(iSD);
iix_tot_sort4 = spikes_samp_new_all{iSD};
spikes_new = iix_tot_sort4/fs*1000;
spikes_new = spikes_new(:).';

if nargout>1
    SNR_avg = nan(size(X,2),size(X,3));
    for j = 1:size(X,3)
        SNR_avg(:,j) = 10*log10( [nanmean(Sm(:,j).^2)*ones(1,size(X,2))]./nanmean([X(:,:,j)-Sm(:,j)*ones(1,size(X,2))].^2) );
    end
end



if do_plot
    rr_samp4 = diff(iix_tot_sort4);
    Nmed = 20;
    if length(rr_samp4) > Nmed*4;
        rr_samp4_med = medfilt1([nanmedian(rr_samp4(1:Nmed))*ones(Nmed/2,1);rr_samp4;nanmedian(rr_samp4(end-Nmed+1:end))*ones(Nmed/2,1)],Nmed);
        rr_samp4_med = rr_samp4_med(Nmed/2+1 : end-Nmed/2);
    else
        rr_samp4_med = nan;
    end
    
    figure,
    ax(1)= subplot(311);
    plot(sig(:,1)),hold on,plot(iix_tot_sort4,sig(iix_tot_sort4,1),'or')
    ax(2)= subplot(312);
    plot(cf3),hold on,plot(iix_tot_sort4,cf3(iix_tot_sort4),'or')
    ax(3)= subplot(313);
    plot(iix_tot_sort4(2:end),diff(iix_tot_sort4),'ok'),
    hold on
    plot(iix_tot_sort4(2:end),rr_samp4_med,'-r')
    plot(spikes_samp(2:end),diff(spikes_samp),'square')
    title(['SDori=',num2str(nanstd(diff(spikes_samp,2)),'%1.2f'),'vs SD=',num2str(nanstd(diff(iix_tot_sort4,2)),'%1.2f')])
    linkaxes(ax,'x')
end




end

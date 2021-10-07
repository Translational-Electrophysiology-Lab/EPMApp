function [RVI,Disp_Local] = RVI_calculation_fun(dt,rt,xyz,Radius_all);

ari = rt-dt;

for ir = 1:length(Radius_all);
    
    rvi_temp.RVImin = nan(size(dt));
    rvi_temp.RVImean = nan(size(dt));
    rvi_temp.RVImean_norm = nan(size(dt));
    rvi_temp.RVImin_norm = nan(size(dt));
    rvi_temp.chan = cell(1,size(dt,2));
    rvi_temp.chan_ok = cell(1,size(dt,2));
    rvi_temp.Nnodes = nan(size(dt));
    rvi_temp.Dist = cell(1,size(dt,2));
    %=
    Disp_Local_temp.range_AT = nan(size(dt));
    Disp_Local_temp.range_RT = nan(size(dt));
    Disp_Local_temp.range_ARI = nan(size(dt));
    
    Disp_Local_temp.grad_max_AT = nan(size(dt));
    Disp_Local_temp.grad_max_RT = nan(size(dt));
    Disp_Local_temp.grad_max_ARI = nan(size(dt));
    
    Disp_Local_temp.grad_mean_AT = nan(size(dt));
    Disp_Local_temp.grad_mean_RT = nan(size(dt));
    Disp_Local_temp.grad_mean_ARI = nan(size(dt));
    
    
    
    Disp_Local_temp.std_AT = nan(size(dt));
    Disp_Local_temp.std_RT = nan(size(dt));
    Disp_Local_temp.std_ARI = nan(size(dt));
    
    Disp_Local_temp.chan = cell(1,size(dt,2));
    Disp_Local_temp.Dist =cell(1,size(dt,2));
    Disp_Local_temp.Nnodes = nan(size(dt));
    
    for ic = 1:size(dt,2)
        d = sqrt(sum( (xyz  - repmat(xyz(ic,:),[size(xyz ,1),1])).^2 ,2));
        [Dim,ii] = sort(d,'ascend');
        ii(1)=[];Dim(1)=[];
%         ii = ii(Dim<=Radius_all(ir));
        ii = ii(Dim<=Radius_all(ir) & Dim>1); % 16/03
        d = d(ii);
        
        if ~isempty(ii)
            RVIall = rt(:,ic)*ones(1,length(ii)) - dt(:,ii);
            %             % RVIi - Delta(AT)
            %             RVIall2 = ari(:,ic)*ones(1,length(ii)) - (dt(:,ii) - dt(:,ic)*ones(1,length(ii)));
            %             % RVIj - Delta(RT)
            %             RVIall3 = ari(:,ii) - (rt(:,ii) - rt(:,ic)*ones(1,length(ii)));
            
            rvi_temp.RVImin(:,ic) =  nanmin(RVIall,[],2);
            rvi_temp.RVImean(:,ic) =  nanmean(RVIall,2);
            
            % normalise by distance
            w = 1./(ones(size(RVIall,1),1)*d');
            iiko = isnan(rt(:,ic) - dt(:,ii));
            w(iiko) = nan;
            rvi_temp.RVImean_norm(:,ic)  =  nanmean((RVIall.*w)./nanmean(w,2),2);
            rvi_temp.RVImin_norm(:,ic)  =  min([(RVIall.*w)./nanmean(w,2)],[],2);
            %
            rvi_temp.chan{ic} = ii;
            rvi_temp.chan_ok{ic} = ~iiko;
            rvi_temp.Dist{ic} = d;
            rvi_temp.Nnodes(:,ic) = sum(~isnan(RVIall),2);
            
            %%
            Disp_Local_temp.range_AT(:,ic) = range(dt(:,[ic;ii]),2);
            Disp_Local_temp.range_RT(:,ic) = range(rt(:,[ic;ii]),2);
            Disp_Local_temp.range_ARI(:,ic) = range(ari(:,[ic;ii]),2);
            
            Disp_Local_temp.std_AT(:,ic) = nanstd(dt(:,[ic;ii]),1,2);
            Disp_Local_temp.std_RT(:,ic) = nanstd(rt(:,[ic;ii]),1,2);
            Disp_Local_temp.std_ARI(:,ic) = nanstd(ari(:,[ic;ii]),1,2);
            
            Disp_Local_temp.grad_max_AT(:,ic) = max(abs(dt(:,ii)-dt(:,ic))./(ones(size(dt,1),1)*d.'),[],2);
            Disp_Local_temp.grad_max_RT(:,ic) = max(abs(rt(:,ii)-rt(:,ic))./(ones(size(dt,1),1)*d.'),[],2);
            Disp_Local_temp.grad_max_ARI(:,ic) = max(abs(ari(:,ii)-ari(:,ic))./(ones(size(dt,1),1)*d.'),[],2);
            
            Disp_Local_temp.grad_mean_AT(:,ic) = nanmean(abs(dt(:,ii)-dt(:,ic))./(ones(size(dt,1),1)*d.'),2);
            Disp_Local_temp.grad_mean_RT(:,ic) = nanmean(abs(rt(:,ii)-rt(:,ic))./(ones(size(dt,1),1)*d.'),2);
            Disp_Local_temp.grad_mean_ARI(:,ic) = nanmean(abs(ari(:,ii)-ari(:,ic))./(ones(size(dt,1),1)*d.'),2);
            
            Disp_Local_temp.chan{ic} = ii;
            Disp_Local_temp.Dist{ic} = d;
            Disp_Local_temp.Nnodes(:,ic) = sum(~isnan(RVIall),2);
        end
    end
    
    name = ['D_',num2str(Radius_all(ir))];name(name=='.')='_';
    RVI.(name)= rvi_temp;
    Disp_Local.(name)= Disp_Local_temp;
end
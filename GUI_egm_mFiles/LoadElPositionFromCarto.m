function [Position,Pos_sens] = LoadElPositionFromCarto(loadname);



fid=fopen(loadname);
if fid==-1
    warning(['File not read: ',loadname])
    
    Position = [];
    return
end
fgetl(fid);fgetl(fid);
dataArray = textscan(fid,'%f %f %f %f %f','CollectOutput',1);
fclose(fid);


if ~isempty(strfind(loadname,'20_POLE_')) % 20 Pole
    
    ii = strfind(loadname,'MAGNETIC_20_POLE');
    name = loadname(ii(1)-7:ii(1)-2);
    ii = strfind(name,'_P');
    Position.point_number = str2double(name(ii(1)+2:end));
    clear ii name
    
    
    ii = [1;find(abs(diff(dataArray{1}(:,1)))>0.5)+1];
    %     if sum(abs(diff(ii,2)))==0
    
    %     chn = [1:length(ii)];
    %     xyz = nan(ii(2)-ii(1),3,length(chn));
    %     t = nan(ii(2)-ii(1),length(chn));
    %     ii = [ii;size(dataArray{1},1)+1];
    %
    %     for i = 1:length(ii)-1
    %         hh = ii(i):ii(i+1)-1;
    %         xyz(hh-hh(1)+1,:,i) =  dataArray{1}(hh,3:5);
    %         t(hh-hh(1)+1,i) =  dataArray{1}(hh,2);
    %     end
    %     D = squeeze(sqrt(sum(xyz.^2,2)));
    
    chn = [1:length(ii)];
    xyz = nan(max(diff(ii)),3,length(chn));
    t = nan(max(diff(ii)),length(chn));
    ii = [ii;size(dataArray{1},1)+1];
    
    if length(ii)<=4
         Position = [];
         Pos_sens = [];
        return;
    end
    
    for i = 1:length(ii)-1
        hh = ii(i):ii(i+1)-1;
        xyz(hh-hh(1)+1,:,i) =  dataArray{1}(hh,3:5);
        t(hh-hh(1)+1,i) =  dataArray{1}(hh,2);
    end
%     D = squeeze(sqrt(sum(xyz.^2,2)));
    Position.t = t(:,3:end);
    Position.xyz = xyz(:,:,3:end); % Exclude sensors (first 2)
%     Position.D = D(:,3:end);
    Position.Labels = [1:size(Position.xyz,3)];
    
    Pos_sens.t = t(:,1:2);
    Pos_sens.xyz = xyz(:,:,1:2);
%     Pos_sens.D = D(:,1:2);
    Pos_sens.Labels = [1:2];
    
    %% PCA
    pca_coeff = nan(3,3,size(xyz,3)); % direction of main components
    pca_v = nan(3,size(xyz,3)); % strength of main components
    
    if size(xyz,1)>1
        for j = 1:size(xyz,3)
            [c,a,v] = pca([xyz(:,:,j)]);
            pca_coeff(:,:,j) = c;
            pca_v(:,j) = v;
        end
    end
    Position.pca_coeff = pca_coeff(:,:,3:end);
    Position.pca_v = pca_v(:,3:end);
    Pos_sens.pca_coeff = pca_coeff(:,:,1:2);
    Pos_sens.pca_v = pca_v(:,1:2);
    
%     % To plot directions on original coordinates
%     j = 1;
%     P1 = [-Position.pca_coeff(1,:,j);Position.pca_coeff(1,:,j)]*sqrt(Position.pca_v(1,j)) + mean([Position.xyz]);
%     P2 = [-Position.pca_coeff(2,:,j);Position.pca_coeff(2,:,j)]*sqrt(Position.pca_v(2,j)) + mean([Position.xyz]);
%     P3 = [-Position.pca_coeff(3,:,j);Position.pca_coeff(3,:,j)]*sqrt(Position.pca_v(3,j)) + mean([Position.xyz]);
%     figure
%     plot3(Position.xyz(:,1,j),Position.xyz(:,2,j),Position.xyz(:,3,j));hold on
%     plot3(P1(:,1),P1(:,2),P1(:,3),'-r');
%     plot3(P2(:,1),P2(:,2),P2(:,3),'color',[1 1 1]*.6);
%     plot3(P3(:,1),P3(:,2),P3(:,3),'color',[1 1 1]*.6);
%     axis(gca,'image')
    
    %% On annotation
    loadname_annot = [loadname(1:end-4),'_OnAnnotation.txt'];
    if exist(loadname_annot,'file');
        fid=fopen(loadname_annot);
        
        fgetl(fid);fgetl(fid);
        dataArray = textscan(fid,'%f %f %f %f %f','CollectOutput',1);
        fclose(fid);
        
        %-
        [~,it] = min([Position.t - nanmedian(dataArray{1}(3:end,2))]);
        Position.t_OnAnnot = it;
        it = nanmedian(it);
        
        Position.xyz_OnAnnot = squeeze(Position.xyz(it,:,:))';
        Pos_sens.xyz_OnAnnot = squeeze(Pos_sens.xyz(it,:,:))';
       
        
        %         Position.t_OnAnnot(1:size(dataArray{1}(1:2,:),1)) = dataArray{1}(1:2,2);
        %         Position.xyz_OnAnnot(1:size(dataArray{1}(1:2,:),1),:) = dataArray{1}(1:2,3:5);
        %         Position.D_OnAnnot(1:size(dataArray{1}(1:2,:),1)) = sqrt(sum(dataArray{1}(1:2,3:5).^2,2));
        %         Position.Labels_OnAnnot(1:size(dataArray{1}(1:2,:),1)) = dataArray{1}(1:2,1);
    end
    
    
elseif ~isempty(strfind(loadname,'NAVISTAR_CONNECTOR_'))% Mapping catheter
    
    ii = strfind(loadname,'NAVISTAR_CONNECTOR_');
    name = loadname(ii(1)-7:ii(1)-2);
    ii = strfind(name,'_P');
    Position.point_number = str2double(name(ii(1)+2:end));
    clear ii name
    
    chn = unique(dataArray{1}(:,1));
    xyz = nan(sum(dataArray{1}(:,1)==1),3,length(chn));
    t = nan(sum(dataArray{1}(:,1)==1),length(chn));
    
    for i = 1:length(chn)
        hh = find(dataArray{1}(:,1)==chn(i));
        xyz(hh-hh(1)+1,:,i) =  dataArray{1}(hh,3:5);
        t(hh-hh(1)+1,i) =  dataArray{1}(hh,2);
    end
    D = squeeze(sqrt(sum(xyz.^2,2)));
    
    iim = min(size(xyz,3),4);
    Position.t = t(:,1:iim);
    Position.xyz = xyz(:,:,1:iim);
%     Position.D = D(:,1:iim);
    Position.Labels = chn(1:iim);
      
    %% PCA
    pca_coeff = nan(3,3,size(xyz,3)); % direction of main components
    pca_v = nan(3,size(xyz,3)); % strength of main components
    for j = 1:size(xyz,3)
        [c,~,v] = pca([xyz(:,:,j)]);
        if ~isempty(c)&&size(c,2)==3
        pca_coeff(:,:,j) = c;
        pca_v(:,j) = v;
        end
    end
    Position.pca_coeff = pca_coeff(:,:,1:iim);
    Position.pca_v = pca_v(:,1:iim);
    
    if length(chn)>4
        Pos_sens.t = t(:,5:size(t,2));
        Pos_sens.xyz = xyz(:,:,5:size(t,2));
        %         Pos_sens.D = D(:,size(t,2));
        Pos_sens.Labels = chn(5:size(t,2));
        Pos_sens.pca_coeff = pca_coeff(:,:,5);
        Pos_sens.pca_v = pca_v(:,5);
        
    else
        Pos_sens.t = [];
        Pos_sens.xyz = [];
        %         Pos_sens.D = [];
        Pos_sens.Labels = [];
        Pos_sens.pca_coeff = [];
        Pos_sens.pca_v = [];
        
    end
    
  
      %% On annotation
    loadname_annot = [loadname(1:end-4),'_OnAnnotation.txt'];
    if exist(loadname_annot,'file');
        fid=fopen(loadname_annot);
        
        fgetl(fid);fgetl(fid);
        dataArray = textscan(fid,'%f %f %f %f %f','CollectOutput',1);
        fclose(fid);
        
        %-
        Position.t_OnAnnot = nan(4,1);
        Position.xyz_OnAnnot = nan(4,3);
%         Position.D_OnAnnot = nan(4,1);
        Position.Labels_OnAnnot = nan(4,1);
        
        Pos_sens.t_OnAnnot = nan(2,1);
        Pos_sens.xyz_OnAnnot = nan(2,3);
%         Pos_sens.D_OnAnnot = nan(2,1);
        Pos_sens.Labels_OnAnnot = nan(2,1);
        
        chn = unique(dataArray{1}(:,1));
        
        Position.t_OnAnnot(dataArray{1}(chn(1:iim),1)) = dataArray{1}(1:iim,2);
        Position.xyz_OnAnnot(dataArray{1}(chn(1:iim),1),:) = dataArray{1}(1:iim,3:5);
%         Position.D_OnAnnot(dataArray{1}(chn(1:iim),1)) = sqrt(sum(dataArray{1}(1:iim,3:5).^2,2));
        Position.Labels_OnAnnot(dataArray{1}(chn(1:iim),1)) = dataArray{1}(1:iim,1);
        
        ll=find(chn>4);
        if ~isempty(ll)
            Pos_sens.t_OnAnnot = dataArray{1}(ll,2);
            Pos_sens.xyz_OnAnnot = dataArray{1}(ll,3:5);
%             Pos_sens.D_OnAnnot = sqrt(sum(dataArray{1}(ll,3:5).^2,2));
            Pos_sens.Labels_OnAnnot = dataArray{1}(ll,1);
        else
            Pos_sens.t_OnAnnot = [];
            Pos_sens.xyz_OnAnnot = [];
%             Pos_sens.D_OnAnnot = [];
            Pos_sens.Labels_OnAnnot = [];
        end

    end
    
else
    error('Connector not found')
end


%% =====================================
%% To check if On_Annotation are correct
% d0 = nan(1,size(Position.xyz,3)); % this should be equal to 1:4 (ablattion catheter) o 1:20 (20A/B) 
% for k = 1:size(Position.xyz,3)
%     for i = 1:size(Position.xyz,3)
%         d(i) = min(sum(abs(Position.xyz(:,:,k)-Position.xyz_OnAnnot(i,:)),2));
%     end
%     j = find(d==0);
%     if ~isempty(j)
%         d0(k)=j;
%     end
% end

%% =====================================
%% To plot catheters + directions


% figure,
% for j = 1 : size(Position.xyz,3)
%     plot3(Position.xyz(:,1,j),Position.xyz(:,2,j),Position.xyz(:,3,j));
%     hold on
%     P1 = [-Position.pca_coeff(1,:,j);Position.pca_coeff(1,:,j)]*sqrt(Position.pca_v(1,j)) + mean(Position.xyz(:,:,j));
%     P2 = [-Position.pca_coeff(2,:,j);Position.pca_coeff(2,:,j)]*sqrt(Position.pca_v(2,j))+ mean(Position.xyz(:,:,j));
%     P3 = [-Position.pca_coeff(3,:,j);Position.pca_coeff(3,:,j)]*sqrt(Position.pca_v(3,j))+ mean(Position.xyz(:,:,j));
%     
%     plot3(P1(:,1),P1(:,2),P1(:,3),'-k');hold on
%     plot3(P2(:,1),P2(:,2),P2(:,3),'color',[1 1 1]*.6);hold on
%     plot3(P3(:,1),P3(:,2),P3(:,3),'color',[1 1 1]*.6);hold on
% end
% axis(gca,'image')








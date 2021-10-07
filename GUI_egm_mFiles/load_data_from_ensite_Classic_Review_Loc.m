function [Virtuals_geo,Geometry,Param] = load_data_from_ensite_Classic_Review_Loc(filename1);

fid=fopen(filename1);
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
Virtuals_geo.Surf_ID = fgetl(fid);
Virtuals_geo.Surf_Name = fgetl(fid);
Virtuals_geo.Surf_Origin = fgetl(fid);
x =fscanf(fid,'%f');
Virtuals_geo.xyz = reshape(x,[size(x,1)/3 3]);
Virtuals_geo.Label = cell(1,size(x,1)/3);
for i = 1:size(x,1)/3
    Virtuals_geo.Label(i) = {num2str(i)};
end

Geometry.Group = fgetl(fid);
Geometry.Group_ID = fgetl(fid);
Geometry.Group_name = fgetl(fid);
a = fgetl(fid);
clear x
x =fscanf(fid,'%f');
Geometry.MESH.vertices = reshape(x,[size(x,1)/3 3]);
a = fgetl(fid);
clear x
x =fscanf(fid,'%f');
Geometry.MESH.faces = reshape(x,[size(x,1)/3 3])+1;

% a = fgetl(fid);
% a = fgetl(fid);
% a = fgetl(fid);
%
% for i =1:256
%     a = fgetl(fid);
%     ii1 = find(a=='(');
%     ii2 =  find(a==')');
%     A{i} = a;
%     xyzEnGuide_1(i,:) =str2num(a(ii1(1)+1:ii2(1)-1));
%     xyzEnGuide_2(i,:) =str2num(a(ii1(2)+1:ii2(2)-1));
% end

Param.Header = Header;
Param.Stduy = Stduy;
Param.Time_start = Time_start;
Param.Time_end = Time_end;
Param.Time_Duration = Time_Duration;



fclose all;


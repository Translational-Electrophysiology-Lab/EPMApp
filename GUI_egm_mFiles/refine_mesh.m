function [tri_new,v_new] = refine_mesh(tri,v);


% 
% tri = geo.Cmesh.triangles;
% v = geo.Cmesh.vertices;

vn = nan(size(tri,1),3);
trin = nan(3*size(tri,1),3);
for i = 1:size(tri,1);
        vn(i,:) = nanmean(v(tri(i,:),:));
        in = size(v,1)+i;
        trin((i-1)*3+1,:) =  [tri(i,1) tri(i,2) in];
        trin((i-1)*3+2,:) =  [tri(i,1) in tri(i,3)];
        trin((i-1)*3+3,:) =  [in tri(i,2) tri(i,3)];
end

tri_new =[trin];
v_new = [v;vn];
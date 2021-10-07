function geo_new = get_Ensite_mesh(xyz);


if size(xyz,1)==256
    N = 256-17;
    clear connectivity
    for i = 1:N
        connectivity(2*i-1,:) = [i i+16 i+17];
        connectivity(2*i,:) = [i i+1 i+17];
    end
    connectivity = [connectivity;1 16 17];
    connectivity = [connectivity;256 256-16 256-15];
    i0=1;
    for i = 1:16-2
        connectivity = [connectivity;i0 i0+i i0+1+i];
    end
    i0=241;
    for i = 1:16-2
        connectivity = [connectivity;i0 i0+i i0+i+1];
    end
    
elseif size(xyz,1)==64
    N = 64-9;
    
    clear connectivity
    
    for i = 1:N
        
        connectivity(2*i-1,:) = [i i+8 i+9];
        
        connectivity(2*i,:) = [i i+1 i+9];
        
    end
    
    connectivity = [connectivity;1 8 9];
    connectivity = [connectivity;64 64-8 64-7];
    
    i0=1;
    for i = 1:8-2
        connectivity = [connectivity;i0 i0+i i0+1+i];
    end
    
    i0=57;
    for i = 1:8-2
        connectivity = [connectivity;i0 i0+i i0+i+1];
    end
    
elseif size(xyz,1)==2048
    N = 2048-65;
    clear connectivity
    for i = 1:N
        connectivity(2*i-1,:) = [i i+64 i+65];
        connectivity(2*i,:) = [i i+1 i+65];
    end
    connectivity = [connectivity;1 64 65];
    connectivity = [connectivity;2048 2048-64 2048-63];
    i0=1;
    for i = 1:64-2
        connectivity = [connectivity;i0 i0+i i0+1+i];
    end
    i0=1985;
    for i = 1:64-2
        connectivity = [connectivity;i0 i0+i i0+i+1];
    end
    
else
    error('Mesh only possible for grid of 64 or 256 nodes')
end


geo_new.Mesh_ventricles = connectivity;
geo_new.Nodes_ventricles =  xyz;

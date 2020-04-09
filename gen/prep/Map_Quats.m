function MapQ=Map_Quats(MapOut)
%generates quaterion maps & symmetric versions
%Quaternion is [w,XYZ]

MapQ.Q_map=zeros(MapOut.ypts,MapOut.xpts,4);
MapQ.Q_map_sym=cell(MapOut.ypts,MapOut.xpts);
for x_co=1:MapOut.xpts
    
    for y_co=1:MapOut.ypts
        
        %extract the crystal type of this pt
        crystal_type=MapOut.Crystal(y_co,x_co);
        
        %convert G to Q
        Q_Pt=conv_G_to_Q(MapOut.GSample(:,:,y_co,x_co));
        MapQ.Q_map(y_co,x_co,:)=Q_Pt;
        
        %calculate symmetries
        MapQ.Q_map_sym{y_co,x_co}=Q_Sym(Q_Pt,crystal_type);
    end
end
MapQ.xpts=MapOut.xpts;
MapQ.ypts=MapOut.ypts;
MapQ.Crystal=MapOut.Crystal;
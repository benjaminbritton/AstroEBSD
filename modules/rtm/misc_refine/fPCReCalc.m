function [MapInfo2] = fPCReCalc(Best_PC,PatternData,MapInfo)
%FPCRECALC Fit a plane to pattern centre measurements

%create the X and Y list from the PNum
num_pts=numel(PatternData.P);

n_r=zeros(num_pts,1);
n_c=zeros(num_pts,1);
n_X_beam=zeros(num_pts,1);
n_Y_beam=zeros(num_pts,1);

for n=1:num_pts
    %find the location of this point & create a list
    [n_r(n),n_c(n)]=find(MapInfo.Data_InputMap.PMap==PatternData.P(n));
    n_X_beam(n)=MapInfo.Data_InputMap.XBeam_Map(n_r(n),n_c(n));
    n_Y_beam(n)=MapInfo.Data_InputMap.YBeam_Map(n_r(n),n_c(n));
end

%fit the data to a plane
RHS=[n_X_beam(:),n_Y_beam(:),0*n_Y_beam(:)+1];

if num_pts > 4
    PC_x_mod=robustfit(RHS,Best_PC(:,1),'bisquare',4.685,0);
    PC_y_mod=robustfit(RHS,Best_PC(:,2),'bisquare',4.685,0);
    PC_z_mod=robustfit(RHS,Best_PC(:,3),'bisquare',4.685,0);
else
    disp('Insufficent number of pts to fit a plane')
end

%use the model to build up the map back
PCX=MapInfo.Data_InputMap.XBeam_Map*PC_x_mod(1)+MapInfo.Data_InputMap.YBeam_Map*PC_x_mod(2)+PC_x_mod(3);
PCY=MapInfo.Data_InputMap.XBeam_Map*PC_y_mod(1)+MapInfo.Data_InputMap.YBeam_Map*PC_y_mod(2)+PC_y_mod(3);
DD= MapInfo.Data_InputMap.XBeam_Map*PC_z_mod(1)+MapInfo.Data_InputMap.YBeam_Map*PC_z_mod(2)+PC_z_mod(3);

%copy over MapInfo
MapInfo2=MapInfo;
MapInfo2.Data_InputMap.PCX=PCX;
MapInfo2.Data_InputMap.PCY=PCY;
MapInfo2.Data_InputMap.DD =DD;
end


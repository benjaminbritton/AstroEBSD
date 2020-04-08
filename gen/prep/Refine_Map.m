%% Fit a better PC model

function [MapRefine]=Refine_Map(no_points,InputUser,Settings_Cor,RTM,file)

close all

%% Load metadata
InputUser.HDF5_file = file;
RTM.Bin_loc=fullfile(InputUser.Astro_loc,'phases','masterpatterns');
RTM.Phase_Folder=fullfile(InputUser.Astro_loc,'phases');

[MapInfo.MapData,MicroscopeData,PhaseData,MapInfo.EBSPData ]=bReadHDF5( InputUser );
[MapInfo.Data_InputMap] = EBSD_Map(MapInfo.MapData,MicroscopeData);

Radon=MapInfo.Data_InputMap.RadonQuality;
SEMImage=transpose(MicroscopeData.SEMImage(:,:,1));

MapSize1=max(MapInfo.MapData.YBeam)-min(MapInfo.MapData.YBeam)+1; % number of rows
MapSize2=max(MapInfo.MapData.XBeam)-min(MapInfo.MapData.XBeam)+1; % number of columns

%% Get user to click some points
figure
I=imagesc(Radon);
colormap('gray')

t=text(8,10,['Select ',num2str(no_points),' points for refinement...'],'Color','blue','FontSize',14);
drawnow

for i = 1:no_points
    p(i) = images.roi.Point;
end

for i = 1:no_points
    draw(p(i))
end

for i =1:no_points
    locations(i,:)=p(1,i).Position;
end


%% Run the PC refinement 
locations=floor(locations);

Refine.ss=0.08; %initial probe volume
Refine.p=2; %order of polynomial to fit to tetrahedron
Refine.n_its=10;
Refine.reindex=1;
Refine.print=0;

for i=1:no_points
   
    delete(t)
    t=text(8,10,['Refining point ',num2str(i),' of ',num2str(no_points),' - do not close this window' ],'Color','red','FontSize',14);
    drawnow
    
    
    PC_start=[MapInfo.Data_InputMap.PCX(locations(i,2),locations(i,1)),MapInfo.Data_InputMap.PCY(locations(i,2),locations(i,1)),MapInfo.Data_InputMap.DD(locations(i,2),locations(i,1))]; %initial value for PC
    Eulers=[MapInfo.Data_InputMap.phi1(locations(i,2),locations(i,1))*degree,MapInfo.Data_InputMap.PHI(locations(i,2),locations(i,1))*degree,MapInfo.Data_InputMap.phi2(locations(i,2),locations(i,1))*degree]; %initial value for Eulers
    
    n=MapInfo.Data_InputMap.PMap(locations(i,2),locations(i,1));
    
    [ RefPat ] = bReadEBSP(MapInfo.EBSPData,n);
    [ RefPat_cor] = EBSP_BGCor( RefPat,Settings_Cor );
    
    [Refine]=PC_refine(Eulers,Refine,PC_start,RefPat_cor,MicroscopeData,RTM,InputUser.Phases{1});
    
    MapRefine.PC_out(i,:)=Refine.PC_out;
    
    MapRefine.PC_in(i,:)=PC_start;
    MapRefine.Eulers_out(i,:)=Refine.Eulers_out;
    MapRefine.Eulers_in(i,:)=Eulers;
    MapRefine.Increase(i)=Refine.Increase;
    MapRefine.Locations=locations;
    
    p(i).Color=[1,0,0];
    
    delete(t)
    t=text(8,10,['Refining point ',num2str(i),' of ',num2str(no_points),' - do not close this window'  ],'Color','red','FontSize',14);
    drawnow
end    
  
%% Fit plane to PCX, PCY, PCZ
PCout=MapRefine.PC_out;
Locs=double(locations)./double([MapSize2,MapSize1]);

%% Want to fit a plane for PCX, PCY, PCZ

%has columns x loc, y loc, PCX
pcx_matrix=[Locs,PCout(:,1)];
pcy_matrix=[Locs,PCout(:,2)];
pcz_matrix=[Locs,PCout(:,3)];

%% fit a plane to PCX
[coeff_x,score,latent]=pca(pcx_matrix,'NumComponents',2,'Centered','on');
v1_x=coeff_x(:,1); %vector 1 in the plane
v2_x=coeff_x(:,2); %vector 2 in the plane
v3_x=cross(v1_x,v2_x); %plane normal

% fit a plane to PCY
[coeff_y,score,latent]=pca(pcy_matrix,'NumComponents',2,'Centered','on');
v1_y=coeff_y(:,1); %vector 1 in the plane
v2_y=coeff_y(:,2); %vector 2 in the plane
v3_y=cross(v1_y,v2_y); %plane normal

% fit a plane to PCZ
[coeff_z,score,latent]=pca(pcz_matrix,'NumComponents',2,'Centered','on');
v1_z=coeff_z(:,1); %vector 1 in the plane
v2_z=coeff_z(:,2); %vector 2 in the plane
v3_z=cross(v1_z,v2_z); %plane normal

%% Generate a map of the new PCX, PCY, PCZ
MapSize1=double(MapSize1);
MapSize2=double(MapSize2);

for j=1:MapSize2
    for i=1:MapSize1

        x=j/MapSize2;
        y=i/MapSize1;

        %point on plane
        point=pcx_matrix(1,:)';

        PCX(i,j)=(dot(point,v3_x)-v3_x(1)*x-v3_x(2)*y)/v3_x(3);
        PCY(i,j)=(dot(point,v3_y)-v3_y(1)*x-v3_y(2)*y)/v3_y(3);
        PCZ(i,j)=(dot(point,v3_z)-v3_z(1)*x-v3_z(2)*y)/v3_z(3);
    end
end

MapRefine.FittedModel.PCX=PCX;
MapRefine.FittedModel.PCY=PCY;
MapRefine.FittedModel.PCZ=PCZ;


close all
end



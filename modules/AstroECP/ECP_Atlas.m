function [fig_info] = ECP_Atlas(Input_Data,ECP_Pat,screen_int,eangs_n,PC_start)
%tool to plot the ECP from a dynamical pattern and show where you are
fig_info=figure;

layout = tiledlayout(1,2);
layout.TileSpacing = 'tight';

ax1 = nexttile([1 1]);


%generate the sterogram
% [ ~,~,~,~,~, RTM_info ] = Phase_Builder_RTM(Input_Data.Phase_Input,Input_Data.Phase_Folder);
% [screen_int] = Cube_Generate(RTM_info.bin_file,RTM_info.isHex);


%% Plot the 2D stereogram
%create the base grid
num_stereo=2001;
stereo_range=1;
stereo_xl=linspace(-stereo_range,stereo_range,num_stereo);
stereo_yl=linspace(-stereo_range,stereo_range,num_stereo);

[stereo_x,stereo_y]=meshgrid(stereo_xl,stereo_yl);

stereo_r=stereo_x.^2+stereo_y.^2;
%normalize to sterogram
stereo_nu=2./(stereo_r+1);

stereo_vx=stereo_nu.*stereo_x;
stereo_vy=stereo_nu.*stereo_y;
stereo_vz=-1+stereo_nu;

stereo_vl=sqrt(stereo_vx.^2+stereo_vy.^2+stereo_vz.^2);
stereo_vx=stereo_vx./stereo_vl;
stereo_vy=stereo_vy./stereo_vl;
stereo_vz=stereo_vz./stereo_vl;

%rotate this sterogram to put the ECP coordinate frame in the middle
eangs_n=eangs_n*pi/180; %convert to degrees
gmatrix=Input_Data.Rz(eangs_n(3))*Input_Data.Rx(eangs_n(2))*Input_Data.Rz(eangs_n(1));
g_dynamics=Input_Data.Rx(pi/2)*Input_Data.Rz(pi/2);
rotmat=g_dynamics*gmatrix;

r_stereo=[stereo_vx(:),stereo_vy(:),stereo_vz(:)];
r_stereo_rot=r_stereo*rotmat';

%sample from the sphere
% [stereo_pat] =
% Cube_Sample(stereo_vx(:),stereo_vy(:),stereo_vz(:),screen_int,0); %in the
% reference frame

[stereo_rot] = Cube_Sample(r_stereo_rot(:,1),r_stereo_rot(:,2),r_stereo_rot(:,3),screen_int,0);

imagesc(stereo_xl,stereo_yl,reshape(stereo_rot,num_stereo,num_stereo));
axis image; axis xy; colormap('gray');

[EBSD_simulation ] = EBSP_Gnom( ECP_Pat,PC_start); %you can change PC_in if you want

% line_top=[EBSD_simulation.xpts_screen(:,1),EBSD_simulation.ypts_screen(:,1)];
% line_bottom=flipud([EBSD_simulation.xpts_screen(:,end),EBSD_simulation.ypts_screen(:,end)]);
% line_right=[EBSD_simulation.xpts_screen(end,:);EBSD_simulation.ypts_screen(end,:)]';
% line_left=flipud([EBSD_simulation.xpts_screen(1,:);EBSD_simulation.ypts_screen(1,:)]');
% %
% line_allx=[line_top(:,1);line_right(:,1);line_bottom(:,1);line_left(:,1)];
% line_ally=[line_top(:,2);line_right(:,2);line_bottom(:,2);line_left(:,2)];
% line_allz=line_ally*0+1;

% hold on;
% plot(line_allx,line_ally,'b','LineWidth',2);
title('Full Stereogram');

%plot the ECP position
%ECP position
hold on;
x_screen_fixm=EBSD_simulation.x_screen(1)*ones(size(EBSD_simulation.y_screen));
x_screen_fixp=EBSD_simulation.x_screen(end)*ones(size(EBSD_simulation.y_screen));
y_screen_fixm=EBSD_simulation.y_screen(1)*ones(size(EBSD_simulation.x_screen));
y_screen_fixp=EBSD_simulation.y_screen(end)*ones(size(EBSD_simulation.x_screen));

% x_screen_perm=[x_screen_fixm,EBSD_simulation.x_screen,x_screen_fixp,fliplr(EBSD_simulation.x_screen)];
% y_screen_perm=[EBSD_simulation.y_screen,y_screen_fixm,fliplr(EBSD_simulation.y_screen),y_screen_fixp];

ecp_1_x=x_screen_fixm;
ecp_1_y=EBSD_simulation.y_screen;

ecp_2_x=EBSD_simulation.x_screen;
ecp_2_y=y_screen_fixm;

ecp_3_x=x_screen_fixp;
ecp_3_y=fliplr(EBSD_simulation.y_screen);

ecp_4_x=fliplr(EBSD_simulation.x_screen);
ecp_4_y=y_screen_fixp;

[ecp_1_xyz]=norm_xy(ecp_1_x,ecp_1_y,rotmat);
[ecp_2_xyz]=norm_xy(ecp_2_x,ecp_2_y,rotmat);
[ecp_3_xyz]=norm_xy(ecp_3_x,ecp_3_y,rotmat);
[ecp_4_xyz]=norm_xy(ecp_4_x,ecp_4_y,rotmat);

[ecp_1_sx,ecp_1_sy]=stereo_norm([ecp_1_x; ecp_1_y; 0*ecp_1_y+1]);
[ecp_2_sx,ecp_2_sy]=stereo_norm([ecp_2_x; ecp_2_y; 0*ecp_2_y+1]);
[ecp_3_sx,ecp_3_sy]=stereo_norm([ecp_3_x; ecp_3_y; 0*ecp_3_y+1]);
[ecp_4_sx,ecp_4_sy]=stereo_norm([ecp_4_x; ecp_4_y; 0*ecp_4_y+1]);



plot(ecp_1_sx,ecp_1_sy,'r','LineWidth',2);
plot(ecp_2_sx,ecp_2_sy,'g','LineWidth',2);
plot(ecp_3_sx,ecp_3_sy,'b','LineWidth',2);
plot(ecp_4_sx,ecp_4_sy,'m','LineWidth',2);

%% Plot the sphere

ax_s2 = nexttile([1 1]);
[sx,sy,sz]=sphere(1000);
[sphere_sim] = Cube_Sample(sx(:),sy(:),sz(:),screen_int,0);

surf(sx,sy,sz,reshape(sphere_sim,size(sx,1),size(sx,2)),'EdgeColor','none');
colormap('gray');
axis off; axis equal;
axis vis3d;
% view(ax_s2,[90 90]); %z+ view
view(ax_s2,[0,1,0])
camroll(ax_s2,-90)
% axis on;
xlabel('x'); ylabel('y'); zlabel('z')
% title('Z+ hemisphere')
hold on;
% plot3(xyz_r_rot(:,1),xyz_r_rot(:,2),xyz_r_rot(:,3),'b');

%overlay the camera/ecp
plot3(ecp_1_xyz(:,1),ecp_1_xyz(:,2),ecp_1_xyz(:,3),'r');
plot3(ecp_2_xyz(:,1),ecp_2_xyz(:,2),ecp_2_xyz(:,3),'g');
plot3(ecp_3_xyz(:,1),ecp_3_xyz(:,2),ecp_3_xyz(:,3),'b');
plot3(ecp_4_xyz(:,1),ecp_4_xyz(:,2),ecp_4_xyz(:,3),'m');

scatter3(1,0,0,10,'r','filled');
scatter3(0,1,0,10,'g','filled');
scatter3(0,0,1,10,'b','filled');
fig_info.WindowState='maximized'; %make this full frame


end

function [xyz_r_rot]=norm_xy(x_screen_perm,y_screen_perm,rotmat)
% a quick function to take something in the gnomonic frame, and turn it
% onto the sphere, and then rotate
z_screen_perm=1+0*y_screen_perm;

xyz_norm=sqrt(x_screen_perm.^2+y_screen_perm.^2+z_screen_perm.^2);
x_screen_n=x_screen_perm./xyz_norm;
y_screen_n=y_screen_perm./xyz_norm;
z_screen_n=z_screen_perm./xyz_norm;


xyz_r=[x_screen_n(:),y_screen_n(:),z_screen_n(:)];
xyz_r_rot=xyz_r*rotmat';
end

function [line_P_rx,line_P_ry]=stereo_norm(line_r2)

[line_r2n]=norm_xy(line_r2(1,:),line_r2(2,:),eye(3));

%convert these points to the stereogram
line_r_lambda=1./(line_r2n(:,3)+1);

%convert from r_exp to the stereogram
line_P_rx=line_r2n(:,1).*line_r_lambda;
line_P_ry=line_r2n(:,2).*line_r_lambda;


end
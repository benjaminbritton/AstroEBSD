%% Index a Single Pattern

%choose a point from a map
x=10; y=10;
Pat_num=MapOut.P_co(y_co,x_co); %row,column, i.e. y,x
PC=[PCFit_out2.PCx_map(y_co,x_co) PCFit_out2.PCy_map(y_co,x_co) PCFit_out2.PCz_map(y_co,x_co)];

%load the pattern
pattern2 = bReadEBSP(EBSPData,Pat_num);
%bg correct
[ EBSP_cor,outputs ] = EBSP_BGCor( pattern2,Settings_Cor );

%generate the geometry
[ EBSD_geometry ] = EBSP_Gnom( detector,PC );

% radon convert & Peak ID
[ Peak_Centre,Peak_Set_ok,Peak_Set,R_EBSP,R_Edge,R_rho,R_theta ] = EBSP_RadHunt( EBSP_cor,Settings_Rad);

% Convert the bands to normal space
[ nhat_gnom] = EBSP_NormConv( Peak_Centre,[outputs.size],PC);

rot_det=eye(3);
%index for all phases
for num_P=1:num_Phases
    [rotdata{num_P},banddata{num_P}]=EBSP_Index(nhat_gnom,Crystal_LUT{num_P},Settings_LUT{num_P}.thresh_trig,Crystal_UCell{num_P},rot_det); %#ok<PFBNS>
end

%%
num_P=1;

figure;

s1=subplot(2,2,1);
imagesc(pattern2); axis equal; colormap('gray'); axis xy; axis equal; axis tight;
title('Raw Input');

s2=subplot(2,2,2);
imagesc(R_theta,R_rho,R_EBSP); colormap('gray'); axis xy; axis equal; axis tight;
hold on; scatter(Peak_Centre(:,1),Peak_Centre(:,2),'r');
xlabel('Theta ( \circ )'); ylabel('Rho');
title('Radon Transform with peaks found');

s3=subplot(2,2,3);
Plot_EBSPNorms( EBSP_cor,EBSD_geometry,nhat_gnom,s3);
title('Background corrected & RT Bands');

s4=subplot(2,2,4);
Plot_EBSPAnnotated( EBSP_cor,EBSD_geometry,nhat_gnom,rotdata{num_P}.detector,Crystal_UCell{num_P},Crystal_LUT{num_P}.family_norm_list,s4);
title('Indexed Pattern');
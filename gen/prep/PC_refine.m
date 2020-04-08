function [Refine]=PC_refine(Eulers,Refine,PC_start,RefPatCor,MicroscopeData,RTM,PhaseInput)

ss=Refine.ss;

%stepsize_x=Refine.ss;
%stepsize_y=Refine.ss;
%stepsize_z=Refine.ss;

p=Refine.p;

phi1=Eulers(1);
PHI=Eulers(2);
phi2=Eulers(3);

xs=[1,3,2,2,2,2,2];
ys=[2,2,1,3,2,2,2];
zs=[2,2,2,2,1,3,2];

p=Refine.p;
n_its=Refine.n_its;
reindex=Refine.reindex;

PC_best=PC_start;
phbest=0;
grad2=[0,0,0];

skipgradcalc=0;

%% Run the iteration
%disp(['Starting refinement'])
for n=1:n_its
phs=zeros(3,3,3);
% Loop over a cube of surrounding PC candidates

if skipgradcalc==0;
for i=1:7

x=xs(i);
y=ys(i);
z=zs(i);
     
x_itval=(x-2).*ss;
y_itval=(y-2).*ss;
z_itval=(z-2).*ss;

PC_refined=[PC_start(1)+x_itval,PC_start(2)+y_itval,PC_start(3)+z_itval];

tilt=MicroscopeData.TotalTilt;
GMat_test=conv_EA_to_G([phi1,PHI,phi2]);
Rx=@(theta)[1 0 0;0 cos(theta) sin(theta);0 -sin(theta) cos(theta)]; %x rotation
Detector_tilt = Rx(tilt);
rotmat=GMat_test*Detector_tilt;

InputPats=RefPatCor;

%% Run the RTM


%phase=1;
%InputUser.Phase_Input=InputUser.Phases(phase);

[ ~,~,~,~,~, RTM_info ] = Phase_Builder_RTM( {PhaseInput},RTM.Phase_Folder);
cs=loadCIF(RTM_info.cif_file);

[ Crystal_UCell,Crystal_Family,Crystal_LUT,Settings_LUT,num_Phases, RTM_info ] = Phase_Builder_RTM( {PhaseInput},RTM.Phase_Folder);
[screen_int,facedata] = Cube_Generate(RTM_info.bin_file,RTM_info.isHex);

[ SettingsXCF, correction, SettingsXCF2 ] = FFT_Filter_settings( RTM.screensize, RTM.LPTsize );

PatternInfo.ScreenWidth=RTM.screensize;
PatternInfo.ScreenHeight=RTM.screensize;
%Set up the screen
[ EBSP_av ] = EBSP_Gnom( PatternInfo,PC_refined );
[ template_pat ] = EBSP_gen( EBSP_av,rotmat,screen_int,0 ); %generate the EBSP for this iteration

%Window, FFT, filter the pattern
%zero mean & fix stdev
template_pat=template_pat-mean(template_pat(:));
template_pat=template_pat./std(template_pat(:));
[template_library,~]  =fROIEx(template_pat,SettingsXCF.hfilter,SettingsXCF.FFTfilter,SettingsXCF.roiloc,SettingsXCF.filters,SettingsXCF.numroi,SettingsXCF.roisize); %FFT the image
       
TemData=zeros(7,1);
PatternCor=InputPats(:,:,1);

Rz=@(theta)[cos(theta) sin(theta) 0;-sin(theta) cos(theta) 0;0 0 1]; %z rotation
Rx=@(theta)[1 0 0;0 cos(theta) sin(theta);0 -sin(theta) cos(theta)]; %x rotation
Ry=@(theta)[cos(theta) 0 sin(theta);0 1 0; -sin(theta) 0 cos(theta)]; %y rotation

%match from library
[rotmat_1st, shift_xFFT, shift_yFFT,rotmat_Peakheight] = index_FFT2(SettingsXCF, PatternCor, rotmat, template_library);
%refine the orientation
[rotmat_best,regout] = refine4(PatternCor,PC_refined,EBSP_av,rotmat_1st,Rx,Ry,Rz,SettingsXCF,SettingsXCF2,correction,screen_int,RTM_info.isHex,RTM.LPTsize,RTM.iterations); 

eang_template=conv_G_to_EA(rotmat_best*inv(Detector_tilt));
TemData(:,1)=[eang_template regout]';

ph=TemData(7,:);
Phi1_Temp=TemData(1,:);
PHI_Temp=TemData(2,:);
Phi2_Temp=TemData(3,:);

%%
phs(x,y,z)=ph;

%save Euler angles of the new trial cube
if x==2 && y==2 && z==2
    phi1_n=Phi1_Temp;
    PHI_n=PHI_Temp;
    phi2_n=Phi2_Temp;
end

end

[maxph,loc]=max(phs(:));

% z_h=mean(reshape(phs([1,3],[1,3],3),2,2),'all');
% z_l=mean(reshape(phs([1,3],[1,3],1),2,2),'all');
% y_h=mean(reshape(phs([1,3],3,[1,3]),2,2),'all');
% y_l=mean(reshape(phs([1,3],1,[1,3]),2,2),'all');
% x_h=mean(reshape(phs(3,[1,3],[1,3]),2,2),'all');
% x_l=mean(reshape(phs(1,[1,3],[1,3]),2,2),'all');

xcol=reshape(phs(:,2,2),1,3);
ycol=reshape(phs(2,:,2),1,3);
zcol=reshape(phs(2,2,:),1,3);

%sum peakheights across the axes and fit a polynomial 
[fitx]=polyfit([-ss,0,ss],xcol,p); %xaxis
[fity]=polyfit([-ss,0,ss],ycol,p); %yaxis
[fitz]=polyfit([-ss,0,ss],zcol,p); %zaxis

%evaluate the gradient at x=0
grad_x=fitx(p);
grad_y=fity(p);
grad_z=fitz(p);

grad=0.1*[grad_x.*ss,grad_y.*ss,grad_z.*ss];%+0.15*grad2;
grad2=grad;

%stepsize_x=10*abs(grad(1));
%stepsize_y=10*abs(grad(2));
%stepsize_z=10*abs(grad(3));

Refine.PC_it(n,:)=PC_start;
Refine.Eulers_it(n,:)=[phi1,PHI,phi2];
Refine.PH_it(n)=phs(2,2,2);

phtrial=phs(2,2,2);

else
    %% do the template match for just one point
    
    PC_refined=PC_start;%[PC_start(1)+x_itval,PC_start(2)+y_itval,PC_start(3)+z_itval];

    tilt=MicroscopeData.TotalTilt;
    GMat_test=conv_EA_to_G([phi1,PHI,phi2]);
    Rx=@(theta)[1 0 0;0 cos(theta) sin(theta);0 -sin(theta) cos(theta)]; %x rotation
    Detector_tilt = Rx(tilt);
    rotmat=GMat_test*Detector_tilt;

    InputPats=RefPatCor;

    %phase=1;
    %InputUser.Phase_Input=InputUser.Phases(phase);
    cs=loadCIF(RTM_info.cif_file);
    [ Crystal_UCell,Crystal_Family,Crystal_LUT,Settings_LUT,num_Phases, RTM_info ] = Phase_Builder_RTM( {PhaseInput},RTM.Phase_Folder);
    [screen_int,facedata] = Cube_Generate(RTM_info.bin_file,RTM_info.isHex);

    [ SettingsXCF, correction, SettingsXCF2 ] = FFT_Filter_settings( RTM.screensize, RTM.LPTsize );

    PatternInfo.ScreenWidth=RTM.screensize;
    PatternInfo.ScreenHeight=RTM.screensize;
    %Set up the screen
    [ EBSP_av ] = EBSP_Gnom( PatternInfo,PC_refined );
    [ template_pat ] = EBSP_gen( EBSP_av,rotmat,screen_int,0 ); %generate the EBSP for this iteration

    %Window, FFT, filter the pattern
    %zero mean & fix stdev
    template_pat=template_pat-mean(template_pat(:));
    template_pat=template_pat./std(template_pat(:));
    [template_library,~]  =fROIEx(template_pat,SettingsXCF.hfilter,SettingsXCF.FFTfilter,SettingsXCF.roiloc,SettingsXCF.filters,SettingsXCF.numroi,SettingsXCF.roisize); %FFT the image

    TemData=zeros(7,1);
    PatternCor=InputPats(:,:,1);

    Rz=@(theta)[cos(theta) sin(theta) 0;-sin(theta) cos(theta) 0;0 0 1]; %z rotation
    Rx=@(theta)[1 0 0;0 cos(theta) sin(theta);0 -sin(theta) cos(theta)]; %x rotation
    Ry=@(theta)[cos(theta) 0 sin(theta);0 1 0; -sin(theta) 0 cos(theta)]; %y rotation

    %match from library
    [rotmat_1st, shift_xFFT, shift_yFFT,rotmat_Peakheight] = index_FFT2(SettingsXCF, PatternCor, rotmat, template_library);
    %refine the orientation
    [rotmat_best,regout] = refine4(PatternCor,PC_refined,EBSP_av,rotmat_1st,Rx,Ry,Rz,SettingsXCF,SettingsXCF2,correction,screen_int,RTM_info.isHex,RTM.LPTsize,RTM.iterations); 

    eang_template=conv_G_to_EA(rotmat_best*inv(Detector_tilt));
    TemData(:,1)=[eang_template regout]';

    ph=TemData(7,:);
    Phi1_Temp=TemData(1,:);
    PHI_Temp=TemData(2,:);
    Phi2_Temp=TemData(3,:);
    
    phtrial=ph;
    phi1_n=Phi1_Temp;
    PHI_n=PHI_Temp;
    phi2_n=Phi2_Temp;  
end

%if there's an improvement:
if n<n_its && phtrial>phbest
    
    %update the current best PC and peak height
    %move the trial PC (for the next iteration)
    PC_best=PC_start;
    phbest=phtrial;
    
    %update Euler angles to central position
    if reindex==1
        phi1=phi1_n;
        PHI=PHI_n;
        phi2=phi2_n;
    end
    
    if skipgradcalc==0
    %set the next PC to follow the gradient
    PC_start=[PC_best(1)+grad(1),PC_best(2)+grad(2),PC_best(3)+grad(3)];
    else
    end
    
    skipgradcalc=0;
    counter=0;

elseif n<n_its %else reduce the step size
    ss=0.9*ss;
    %and move slightly randomly
    reduction=(5-mod(counter,5))/5;
    
    grad=reduction.*0.01.*[(rand-0.5),(rand-0.5),(rand-0.5)].*PC_start; % -0.25 to 0.25 percent of PC_start
    PC_start=[PC_best(1)+grad(1),PC_best(2)+grad(2),PC_best(3)+grad(3)];
    skipgradcalc=1;
    
    counter=counter+1;
    
end

if Refine.print==1
    disp(['Completed iteration ',num2str(n),': Trial PH: ',num2str(phtrial),', Best PH: ',num2str(phbest)]);%,', next SS: ',num2str(stepsize)])
end
%disp(grad);
%disp(PC_start)

end

[Refine.PH_out,loc]=(max(Refine.PH_it));
Refine.PC_out=Refine.PC_it(loc,:);
Refine.Eulers_out=Refine.Eulers_it(loc,:);
Refine.Pattern=RefPatCor;

Refine.Increase=(Refine.PH_out-Refine.PH_it(1))./Refine.PH_it(1);

end
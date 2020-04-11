function [rotmat_best,RegFinal] = refine5(EBSP_ref,EBSP_geometry,PC_average,rotmat_in,SettingsXCF,screen_int,isHex,RTM_setup)
%[ corrCo, best_rotmat ] = refine2( Detector, i, screensize, R_x, R_y, rotmat_best_SO3, shift_x, shift_y, XCF_type, f, EBSP_ref, LPTsize, correction, SettingsXCF2, Rz, SettingsXCF, Detector2, iterations )
%REFINE Refine the orientation
screensize=EBSP_geometry.size(1);
iterations=RTM_setup.iterations;
LPTsize=RTM_setup.LPTsize;

PC_pat=EBSP_geometry.PC;
%PCs are w.r.t. square patterns

%Work out shift in x and y and also in pixels
shift = PC_pat - PC_average; %actual - average
x_shift = shift(1);
y_shift = shift(2);

%calculate the angle from the X and Y PC shifts
x_angle = atan(y_shift/PC_pat(3)); %x angle uses y shift as a shift in y is equivelent to a rotation in x
y_angle = atan(x_shift/PC_pat(3)); %y angle uses x shift as a shift in x is equivelelent to a rotation in y

%calculate the rotation matricies
rotmat_PCX_correction = RTM_setup.Rx(x_angle);
rotmat_PCY_correction = RTM_setup.Ry(y_angle);

%update the rotation matrix for the pattern
rotmat_best_SO3 = rotmat_PCX_correction * rotmat_PCY_correction * rotmat_in;

FFTData_Ref=EBSP_ref.FFT;

% %FFT the experimental EBSP
% [FFTData_Ref,~]  =fROIEx2(EBSP_cor,SettingsXCF);

%start the iterations
rotmat_best=rotmat_best_SO3;

%LPT settings
rmin=10;
rmax=screensize/sqrt(2);

%unpack settings XCF2

SettingsXCF2.roisize=SettingsXCF.roisize2;
SettingsXCF2.roiloc=SettingsXCF.roiloc2;
SettingsXCF2.numroi=SettingsXCF.numroi2;
SettingsXCF2.filters=SettingsXCF.filters2;
SettingsXCF2.FFTfilter=SettingsXCF.FFtfilter2;
SettingsXCF2.hfilter=SettingsXCF.hfilter2;
SettingsXCF2.mesh=SettingsXCF.mesh;

correction.xlim=SettingsXCF.xlim;
correction.xrange=SettingsXCF.xrange;
correction.yroi_num=SettingsXCF.yroi_num;
correction.ycen=SettingsXCF.ycen;

p=1;

while p < iterations
    
    %cross correlate this pattern in X and Y to get shifts
    %generate a new pattern
    [ EBSP_iter ] = EBSP_gen( EBSP_geometry,rotmat_best,screen_int,isHex ); %generate the EBSP for this iteration
    
    %FFT it
    %     [EBSP_FFT,XCF_data_fill]  =fROIEx2(EBSP_iter,SettingsXCF);
    [EBSP_iter_r,XCF_data_fill] = refine_prep(EBSP_iter,SettingsXCF,RTM_setup);
    
    %     [FFTData_iter,XCF_data_fill]  =fROIEx2(EBSP_iter,SettingsXCF);
    
    %Cross correlate to get an x and y shift
    RegOut1 = fReg( FFTData_Ref,EBSP_iter_r.FFT,SettingsXCF.roisize,SettingsXCF.mesh,XCF_data_fill); %RegOut = [Xshift, Yshift, fullXCFheight, normXCFheight]
    shift_x_new = RegOut1(1);
    shift_y_new = RegOut1(2);
    
    %calculate the new angular update
    angle.xcorrFFT  = atan(shift_y_new/(PC_pat(3) * EBSP_geometry.size(2)));
    angle.ycorrFFT  = atan(shift_x_new/(PC_pat(3) * EBSP_geometry.size(1)));
    
    %Put them into matrix form
    rotmat.xcorrFFT(:,:) = RTM_setup.Rx(-angle.xcorrFFT );
    rotmat.ycorrFFT(:,:) = RTM_setup.Ry(-angle.ycorrFFT );
    
    %join these two together and update best
    rotmat_best=rotmat_best * rotmat.ycorrFFT*rotmat.xcorrFFT;
    
    %log-polar
    %     EBSP_logp = logsample(EBSP_halfcorr, rmin, rmax, screensize/2, screensize/2, RTM_setup.LPTsize, RTM_setup.LPTsize); %Transform the reference image into LPT space, logsample is in logsample
    %Work out z rotation
    % %           perform LPT FFT XCF to get z rotation
    %     EBSP_logref = logsample(EBSP_cor, rmin, rmax, screensize/2, screensize/2, LPTsize, LPTsize); %Transform the reference image into LPT space, logsample is in logsample
    
    %generate the template again and match
    [ EBSP_halfcorr ] = EBSP_gen(EBSP_geometry,rotmat_best,screen_int,isHex );
    EBSP_loghalfcorr = logsample(EBSP_halfcorr, rmin, rmax, screensize/2, screensize/2, LPTsize, LPTsize); %Thransform the x/y corrected image into LPT space
    %measure the Z rotation
    [ rotmat_z_correction] = Angle_ShiftZ( EBSP_ref.logp,EBSP_loghalfcorr,correction,LPTsize,SettingsXCF2,RTM_setup.Rz ); %This generate the Z correction roation matrix for the x/y corrected image
    rotmat_best = rotmat_best * rotmat_z_correction;
    
    %do this without generation of the template again
%     [ rotmat_z_correction] = Angle_ShiftZ( EBSP_ref.logp,EBSP_iter_r.logp,correction,LPTsize,SettingsXCF2,RTM_setup.Rz ); %This generate the Z correction roation matrix for the x/y corrected image
%     rotmat_best = rotmat_best * rotmat.ycorrFFT*rotmat.xcorrFFT * rotmat_z_correction;
    angle.z_up=acos(rotmat_z_correction(1));
    angle_tot=sum(abs([angle.z_up,angle.xcorrFFT,angle.ycorrFFT]));
    
    %close the loop early if needed
    if angle_tot < 8.7266e-04 % if the changes are smaller than 0.05*180/pi , close the iterations
        p=iterations;
    else
        p=p+1;
    end
end

%generate a new pattern
[ EBSP_final ] = EBSP_gen( EBSP_geometry,rotmat_best,screen_int,isHex ); %generate the EBSP for this iteration

%FFT it
[FFTData_final]  =fROIEx(EBSP_final,SettingsXCF.hfilter,SettingsXCF.FFTfilter,SettingsXCF.roiloc,SettingsXCF.filters,SettingsXCF.numroi,SettingsXCF.roisize);

%Cross correlate to get an x and y shift
RegFinal = fReg( FFTData_Ref,FFTData_final,SettingsXCF.roisize,SettingsXCF.mesh,XCF_data_fill); %RegOut = [Xshift, Yshift, fullXCFheight, normXCFheight]

end


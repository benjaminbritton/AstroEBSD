function [rotmat_best,RegFinal] = refine4(EBSP_cor,PC_av,EBSP,rotmat_best_SO3_in,Rx,Ry,Rz,SettingsXCF,SettingsXCF2,correction,screen_int,isHex,LPTsize,iterations)
%[ corrCo, best_rotmat ] = refine2( Detector, i, screensize, R_x, R_y, rotmat_best_SO3, shift_x, shift_y, XCF_type, f, EBSP_ref, LPTsize, correction, SettingsXCF2, Rz, SettingsXCF, Detector2, iterations )
%REFINE Refine the orientation
screensize=EBSP.size(1);

PC_pat=EBSP.PC;
%PCs are w.r.t. square patterns

%Work out shift in x and y and also in pixels
shift = PC_pat - PC_av; %actual - average
x_shift = shift(1);
y_shift = shift(2);

%calculate the angle from the X and Y PC shifts
x_angle = atan(y_shift/PC_pat(3)); %x angle uses y shift as a shift in y is equivelent to a rotation in x
y_angle = atan(x_shift/PC_pat(3)); %y angle uses x shift as a shift in x is equivelelent to a rotation in y

%calculate the rotation matricies
rotmat_PCX_correction = Rx(x_angle);
rotmat_PCY_correction = Ry(y_angle);

%update the rotation matrix for the pattern
rotmat_best_SO3 = rotmat_PCX_correction * rotmat_PCY_correction * rotmat_best_SO3_in;

%FFT the experimental EBSP
[FFTData_Ref,~]  =fROIEx2(EBSP_cor,SettingsXCF);

%start the iterations
rotmat_best=rotmat_best_SO3;

%LPT settings
rmin=10;
rmax=screensize/sqrt(2);

for p= 1:iterations
    %cross correlate this pattern in X and Y to get shifts
    %generate a new pattern
    [ EBSP_iter ] = EBSP_gen( EBSP,rotmat_best,screen_int,isHex ); %generate the EBSP for this iteration
    
    %FFT it
    [FFTData_iter,XCF_data_fill]  =fROIEx2(EBSP_iter,SettingsXCF);
    
    %Cross correlate to get an x and y shift
    RegOut1 = fReg( FFTData_Ref,FFTData_iter,SettingsXCF.roisize,SettingsXCF.mesh,XCF_data_fill); %RegOut = [Xshift, Yshift, fullXCFheight, normXCFheight]
    shift_x_new = RegOut1(1);
    shift_y_new = RegOut1(2);
    
    %calculate the new angular update
    angle.xcorrFFT  = atan(shift_y_new/(PC_pat(3) * EBSP.size(2)));
    angle.ycorrFFT  = atan(shift_x_new/(PC_pat(3) * EBSP.size(1)));
    
    %Put them into matrix form
    rotmat.xcorrFFT(:,:) = Rx(-angle.xcorrFFT );
    rotmat.ycorrFFT(:,:) = Ry(-angle.ycorrFFT );
    
    %join these two together
    rotmat_halfcorr = rotmat.ycorrFFT*rotmat.xcorrFFT;
    rotmat_best=rotmat_best * rotmat_halfcorr;
    
    %update the pattern
    [ EBSP_halfcorr ] = EBSP_gen(EBSP,rotmat_best,screen_int,isHex );
    
    %Work out z rotation
    % %           perform LPT FFT XCF to get z rotation
    EBSP_logref = logsample(EBSP_cor, rmin, rmax, screensize/2, screensize/2, LPTsize, LPTsize); %Transform the reference image into LPT space, logsample is in logsample
    EBSP_loghalfcorr = logsample(EBSP_halfcorr, rmin, rmax, screensize/2, screensize/2, LPTsize, LPTsize); %Thransform the x/y corrected image into LPT space
    
    %measure the Z rotation
    [ rotmat_z_correction] = Angle_ShiftZ( EBSP_logref,EBSP_loghalfcorr,correction,LPTsize,SettingsXCF2,Rz ); %This generate the Z correction roation matrix for the x/y corrected image
    
    %Create a correction rotation matrix
    %             rotmat_full_correction =  rotmat_z_correction * rotmat_halfcorr;  %combine the x, y and z correction matrices to apply to the SO(3) search
    rotmat_best = rotmat_best * rotmat_z_correction; %combine the full correction with the SO(3) search rotation matrix
    
end

%generate a new pattern
[ EBSP_final ] = EBSP_gen( EBSP,rotmat_best,screen_int,isHex ); %generate the EBSP for this iteration

%FFT it
[FFTData_final]  =fROIEx(EBSP_final,SettingsXCF.hfilter,SettingsXCF.FFTfilter,SettingsXCF.roiloc,SettingsXCF.filters,SettingsXCF.numroi,SettingsXCF.roisize);

%Cross correlate to get an x and y shift
RegFinal = fReg( FFTData_Ref,FFTData_final,SettingsXCF.roisize,SettingsXCF.mesh,XCF_data_fill); %RegOut = [Xshift, Yshift, fullXCFheight, normXCFheight]

end


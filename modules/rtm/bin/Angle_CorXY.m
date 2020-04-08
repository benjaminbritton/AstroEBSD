function [ rotmat_XY ] = Angle_CorXY( Detector,shift_xFFT,shift_yFFT,R_x,R_y, EBSP_ref, SettingsXCF, rotmat_bestNDP, f , XCF_type, screensize, i )
%ANGLE_CORXY Formulate the XY rotation corrections based upon the X and Y
%shifts from the XCF

% This code is copyright Alex Foden and Ben Britton 09/04/2019
%
% This RTM-EBSD (refined template indexing) code and its associated scripts
% may only be shared with express and direct permission of
% Alex Foden & Ben Britton
% 
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
% OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
% NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
% BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN 
% ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
% CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

% Requirements:
% MATLAB R2018a or above
% MTEX version 5.2.beta2 or above
% Created by Alex Foden and Ben Britton 28/03/2019
% If you are using a CIF file not in the MTEX toolbox, you will need to add
% the full file path to the cif file to the phase file you are using

switch XCF_type
    case 1
        SettingsXCF.mesh = 3;
        %FFT the refernece image
        [FFTData_Ref,XCF_data_fill]  =fROIEx(EBSP_ref,SettingsXCF.hfilter,SettingsXCF.FFTfilter,SettingsXCF.roiloc,SettingsXCF.filters,SettingsXCF.numroi,SettingsXCF.roisize);
        
        EBSP_bestNDP = image_gen2(Detector, rotmat_bestNDP, f, i);
        %FFT the test image
        [FFTData_test,~]  =fROIEx(EBSP_bestNDP,SettingsXCF.hfilter,SettingsXCF.FFTfilter,SettingsXCF.roiloc,SettingsXCF.filters,SettingsXCF.numroi,SettingsXCF.roisize);
        
        %Cross correlate to get an x and y shift
        RegOut1 = fReg( FFTData_Ref,FFTData_test,SettingsXCF.roisize,SettingsXCF.mesh,XCF_data_fill); %RegOut = [Xshift, Yshift, fullXCFheight, normXCFheight]
        
        shift_xNDP = RegOut1(1);
        shift_yNDP = RegOut1(2);
        
        %Work out the x and y correction rotations
        
        angle.xcorrNDP = atan(shift_yNDP/(Detector.DD(i) * screensize));
        angle.ycorrNDP = atan(shift_xNDP/(Detector.DD(i) * screensize));
        
        %Put them into matrix form
        
        rotmat_xcorrNDP = R_x(angle.xcorrNDP);
        rotmat_ycorrNDP = R_y(angle.ycorrNDP);
        
        rotmat_XY = rotmat_ycorrNDP*rotmat_xcorrNDP;
        
        
    case 2
        angle.xcorrFFT  = atan(shift_yFFT/(Detector.DD(i) * Detector.ScreenWidth));
        angle.ycorrFFT  = atan(shift_xFFT/(Detector.DD(i) * Detector.ScreenWidth));
        
        %Put them into matrix form
        
        rotmat.xcorrFFT(:,:) = R_x(angle.xcorrFFT );
        rotmat.ycorrFFT(:,:) =  R_y(angle.ycorrFFT );
        
        rotmat_XY = rotmat.ycorrFFT* rotmat.xcorrFFT;
        
end
end


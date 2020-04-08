function [ Test_image ] = Library_Gen(EBSP,screen_int,isHex,G_all,XCF_type,SettingsXCF)
% Detector,G_all,SettingsXCF,f, XCF_type )
%LIBRARY_GEN Generates a library of patterns, and FFTs them

% This code is copyright Alex Foden and Ben Britton 09/04/2019
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

%build the detector positions
r = [EBSP.xpts_screen(:), EBSP.ypts_screen(:), EBSP.ypts_screen(:)*0+1].*1./sqrt((EBSP.xpts_screen(:).^2+EBSP.ypts_screen(:).^2+1));
sy=EBSP.size(1);
sx=EBSP.size(2);

switch XCF_type
    case 1 %NDP library (realspace)
        Test_image=zeros(length(EBSP.y_screen),length(EBSP.x_screen),size(G_all,3));
        parfor n = 1:size(G_all,3)
            %rotate the screen
            r2 = r*G_all(:,:,n)';
            %sample the pattern from the interpolant
            [i_data] = Cube_Sample(r2(:,1),r2(:,2),r2(:,3),screen_int,isHex);
            %reshape the output
            Test_image_temp=reshape(i_data,sy,sx);
            
            Test_image_temp=Test_image_temp-mean(Test_image_temp(:));
            Test_image(:,:,n)=Test_image_temp./std(Test_image_temp(:));
        end
        
    case 2 %FFT library (Fourier space)
        %Make and FFT all images on SO(3)
        
        % INCLUDE A PREALLOCATION STEP
        
        parfor n=1:size(G_all,3)
            %rotate the screen
            r2 = r*G_all(:,:,n)';
            %sample the pattern from the interpolant
            [i_data] = Cube_Sample(r2(:,1),r2(:,2),r2(:,3),screen_int,isHex);
            %reshape the output
            P_test_temp=reshape(i_data,sy,sx);
            
            P_test_temp=P_test_temp-mean(P_test_temp(:));
            P_test=P_test_temp./std(P_test_temp(:));
            
            %Window, FFT, filter the patterns
            [Test_image(:,:,n),~]  =fROIEx2(P_test,SettingsXCF); %FFT the image
        end
end

end


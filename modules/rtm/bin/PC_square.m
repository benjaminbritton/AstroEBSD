function [ Data_InputMap2 ] = PC_square( EBSD_DataInfo, Data_InputMap,Settings_Cor )

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

Data_InputMap2=Data_InputMap;

if ~isfield(Settings_Cor,'SquareCrop')
    Settings_Cor.SquareCrop=0;
end

if Settings_Cor.SquareCrop == 1
    %PC_CROPPING_CORRECTION Summary of this function goes here
    %   Detailed explanation goes here
    
    %Add in variables needed for the PC update
    Raw_x_Screensize = EBSD_DataInfo.PW;
    Raw_y_screensize = EBSD_DataInfo.PH;
    
        
    for xx=1:Data_InputMap.xpts
        for yy=1:Data_InputMap.ypts
            Original_Pattern_Centre_x = Data_InputMap.PCX(yy,xx);
            
            Data_InputMap2.PCX(yy,xx)= (Original_Pattern_Centre_x*Raw_x_Screensize-(Raw_x_Screensize-Raw_y_screensize)/2)/(Raw_y_screensize);
            
        end
    end



end

end


% This code is copyright Alex Foden and Ben Britton 09/04/2019
% Do not distribute.
%
% This RTM-EBSD (refined template indexing) code and its associated scripts
% may only be shared with express and direct permission of
% Alex Foden & Ben Britton
% 
% If you would like a copy or access to the repository
% Please contact b.britton@imperial.ac.uk or a.foden16@imperial.ac.uk
%
% Ben Britton has a list of authorised users of this code
% 
% Do not fork
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

%% Write data to a h5 file

%create the file container
OutputUser=InputUser;

% OutputUser.HDF5_file=[OutputUser.EBSD_File(1:end-3) '_RTM.h5']; %If input
% file ends with ".h5"
OutputUser.HDF5_file=[OutputUser.EBSD_File '_Refined_hough.h5']; %If input file doesn't end with ".h5"
%we assume that this is an h5 file, hdf5 will cause issues with the (end-3)
% OutputUser.DataName=InputUser.EBSD_File(1:end-3);
OutputUser.DataName=InputUser.EBSD_File;
OutputUser.HDF5FullFile=fullfile(OutputUser.HDF5_folder,OutputUser.HDF5_file);

%export the data to be a h5 file that looks like a Bruker file
dtype='/EBSD/Data/'; %EBSD data location
htype='/EBSD/Header/'; %Header data location (e.g. microscope settings)

%Pattern centre
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,dtype,'DD',Data_InputMap.DD(:)); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,dtype,'PCX',Data_InputMap.PCX(:)); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,dtype,'PCY',Data_InputMap.PCY(:)); %write the data

%construct arrays from the indexed data
RTM_Phi1=zeros(MicroscopeData.NPoints,1);
RTM_PHI=zeros(MicroscopeData.NPoints,1);
RTM_Phi2=zeros(MicroscopeData.NPoints,1);
RTM_XBeam=RTM_Phi1;
RTM_YBeam=RTM_Phi1;
RTM_XSample=RTM_Phi1;
RTM_YSample=RTM_Phi1;



for p=1:MicroscopeData.NPoints
    RTM_Phi1(p)=TemData(1,p)*180/pi; %save in degrees
    RTM_PHI(p)=TemData(2,p)*180/pi;
    RTM_Phi2(p)=TemData(3,p)*180/pi;
    
    RTM_peakheight(p)=TemData(7,p);
    
    RTM_XBeam(p)=double(Data_InputMap_Start.XBeam_Map(p))-1;
    RTM_YBeam(p)=double(Data_InputMap_Start.YBeam_Map(p))-1;
    RTM_XSample(p)=double(Data_InputMap_Start.XSample(p))-1;
    RTM_YSample(p)=double(Data_InputMap_Start.YSample(p))-1;
end

RTM_peakheight=RTM_peakheight(:); %form as a column 
%orientation data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,dtype,'PHI',RTM_PHI); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,dtype,'phi2',RTM_Phi2); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,dtype,'phi1',RTM_Phi1); %write the data

%beam data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,dtype,'X BEAM',RTM_XBeam); %subtract 1 from the data because of the 0 indexing
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,dtype,'Y BEAM',RTM_YBeam); %subtract 1 from the data because of the 0 indexing
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,dtype,'X SAMPLE',RTM_XSample); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,dtype,'Y SAMPLE',RTM_YSample); %write the data

%indexing data

h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,dtype,'MAD',Data_InputMap.MAD); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,dtype,'MADPhase',Data_InputMap.MADPhase); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,dtype,'NIndexedBands',Data_InputMap.NIndexedBands); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,dtype,'RadonQuality',Data_InputMap.RadonQuality); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,dtype,'PeakHeight',RTM_peakheight); %write the data

%header
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'CameraTilt',MicroscopeData.CameraTilt); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'SampleTilt',MicroscopeData.SampleTilt); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'TotalTilt',MicroscopeData.TotalTilt); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'KV',MicroscopeData.KV); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'NCOLS',MicroscopeData.NCOLS); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'NROWS',MicroscopeData.NROWS); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'XSTEP',MicroscopeData.XSTEP); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'YSTEP',MicroscopeData.YSTEP); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'NPoints',MicroscopeData.NPoints); %write the data

h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'Magnification',MicroscopeData.Magnification); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'WD',MicroscopeData.WD); %write the data

pTime(['Saved to ' OutputUser.HDF5FullFile],time1);
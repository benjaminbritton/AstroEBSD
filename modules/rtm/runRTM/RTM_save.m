function[OutputUser] = RTM_save(RTM,ntiles,InputUser,MapInfo,MicroscopeData)

% This code is copyright Alex Foden and Ben Britton 09/04/2019
% Do not distribute.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
% OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
% NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
% BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN 
% ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
% CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

% Modified by TPM 20/02/2020 for AstroEBSD v2

%% Write data to a h5 file
   

    phase=RTM.Output.phase;
    RTM_Phi1_in=RTM.Output.euler1*180/pi;
    RTM_PHI_in=RTM.Output.euler2*180/pi;
    RTM_Phi2_in=RTM.Output.euler3*180/pi;
    RTM_peakheight_in=RTM.Output.qual;
    
    grid=sqrt(ntiles);


%create the file container
OutputUser=InputUser;

OutputUser.HDF5_file=[OutputUser.EBSD_File(1:end-3) '_RTM.h5']; %If input
% file ends with ".h5"
%OutputUser.HDF5_file=[OutputUser.EBSD_File '_RTM.h5']; %If input file doesn't end with ".h5"
%we assume that this is an h5 file, hdf5 will cause issues with the (end-3)
OutputUser.DataName=InputUser.EBSD_File(1:end-3);
%OutputUser.DataName=InputUser.EBSD_File;

if isfield(OutputUser,'ResultsDir')==0
    OutputUser.ResultsDir=OutputUser.HDF5_folder;
else
    OutputUser.HDF5_folder=OutputUser.ResultsDir;
end

OutputUser.HDF5FullFile=fullfile(OutputUser.HDF5_folder,OutputUser.HDF5_file);

%Delete an HDF5 file of the same name
if exist(OutputUser.HDF5FullFile)==2
    delete(OutputUser.HDF5FullFile)
end

%export the data to be a h5 file that looks like a Bruker file
dtype='/EBSD/Data/'; %EBSD data location
htype='/EBSD/Header/'; %Header data location (e.g. microscope settings)

DD=MapInfo.Data_InputMap.DD(1:grid*MapInfo.cropped_height,1:grid*MapInfo.cropped_width);
PCX=MapInfo.Data_InputMap.PCX(1:grid*MapInfo.cropped_height,1:grid*MapInfo.cropped_width);
PCY=MapInfo.Data_InputMap.PCY(1:grid*MapInfo.cropped_height,1:grid*MapInfo.cropped_width);
XBeam_Map=MapInfo.Data_InputMap.XBeam_Map(1:grid*MapInfo.cropped_height,1:grid*MapInfo.cropped_width);
YBeam_Map=MapInfo.Data_InputMap.YBeam_Map(1:grid*MapInfo.cropped_height,1:grid*MapInfo.cropped_width);
XSample=MapInfo.Data_InputMap.XSample(1:grid*MapInfo.cropped_height,1:grid*MapInfo.cropped_width);
YSample=MapInfo.Data_InputMap.YSample(1:grid*MapInfo.cropped_height,1:grid*MapInfo.cropped_width);

%Pattern centre
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,dtype,'DD',DD(:)); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,dtype,'PCX',PCX(:)); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,dtype,'PCY',PCY(:)); %write the data

%construct arrays from the indexed data
RTM_Phi1=zeros(ntiles*MapInfo.cropped_width*MapInfo.cropped_height,1);
RTM_PHI=zeros(ntiles*MapInfo.cropped_width*MapInfo.cropped_height,1);
RTM_Phi2=zeros(ntiles*MapInfo.cropped_width*MapInfo.cropped_height,1);
RTM_XBeam=RTM_Phi1;
RTM_YBeam=RTM_Phi1;
RTM_XSample=RTM_Phi1;
RTM_YSample=RTM_Phi1;
RTM_peakheight=zeros(ntiles*MapInfo.cropped_width*MapInfo.cropped_height,1);

%RTM_peakheight=RTM_peakheight(:); %form as a column 

for p=1:ntiles*MapInfo.cropped_width*MapInfo.cropped_height
    RTM_Phi1(p)=RTM_Phi1_in(p); %save in degrees
    RTM_PHI(p)=RTM_PHI_in(p);
    RTM_Phi2(p)=RTM_Phi2_in(p);
    
    RTM_peakheight(p)=RTM_peakheight_in(p);
    
    RTM_XBeam(p)=double(XBeam_Map(p))-1;
    RTM_YBeam(p)=double(YBeam_Map(p))-1;
    RTM_XSample(p)=double(XSample(p))-1;
    RTM_YSample(p)=double(YSample(p))-1;
    %disp(p)
end

%RTM_peakheight=RTM_peakheight';

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

h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,dtype,'MADPhase',RTM_Phi1*0+1); %write the data

h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,dtype,'PeakHeight',RTM_peakheight); %write the data

%header
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'CameraTilt',MicroscopeData.CameraTilt); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'SampleTilt',MicroscopeData.SampleTilt); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'TotalTilt',MicroscopeData.TotalTilt); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'KV',MicroscopeData.KV); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'NCOLS',MapInfo.cropped_width*ntiles); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'NROWS',MapInfo.cropped_height*ntiles); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'XSTEP',MicroscopeData.XSTEP); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'YSTEP',MicroscopeData.YSTEP); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'NPoints',ntiles*MapInfo.cropped_width*MapInfo.cropped_height); %write the data

h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'Magnification',MicroscopeData.Magnification); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'WD',MicroscopeData.WD); %write the data
end
%pTime(['Saved to ' OutputUser.HDF5FullFile],time1);
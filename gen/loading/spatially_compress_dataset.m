
%%
InputUser.HDF5_folder='I:\TomMcA\CarbideSet'; %folder where the h5/bcf is located
InputUser.HDF5_file='LR12ht_8kX_1.h5';

[MapData,MicroscopeData,PhaseData,EBSPData ]=bReadHDF5( InputUser );
EBSPData_original=EBSPData
[DataInputMap] = EBSD_Map(MapData,MicroscopeData);

factor=3;           
channum=2048;
%function [DataInputMap2]=spatial_downsample(Data_InputMap,factor);

%%
col_locs=[1:factor:DataInputMap.xpts];
row_locs=[1:factor:DataInputMap.ypts];
[new_cols,new_rows]=meshgrid(col_locs,row_locs);
newrows=reshape(new_rows,1,[]);
newcols=reshape(new_cols,1,[]);

linears=sub2ind([DataInputMap.ypts,DataInputMap.xpts],newrows,newcols);
[new_ypts,new_xpts]=size(new_rows);

%%
fields={'XBeam_Map','YBeam_Map','DD','MAD','MADPhase','NIndexedBands','PCX','PCY','PHI','phi2','phi1','RadonQuality','XSample','YSample'};
for i=1:length(fields)
    field=fields{i};
    %DataInputMap2.(field)=reshape(DataInputMap.(field)(linears),new_ypts,new_xpts);
    DataInputMap2.(field)=reshape(DataInputMap.(field)(linears),new_ypts,new_xpts);
    DataInputMap2.(field)=reshape(DataInputMap2.(field)',[],1);
end

DataInputMap2.XStep=factor*DataInputMap.XStep;
DataInputMap2.YStep=factor*DataInputMap.YStep;
DataInputMap2.X_bline=[1:1:new_xpts];
DataInputMap2.Y_bline=[1:1:new_ypts];
DataInputMap2.X_axis=DataInputMap.X_axis(DataInputMap.X_bline(col_locs));
DataInputMap2.Y_axis=DataInputMap.Y_axis(DataInputMap.Y_bline(row_locs));
DataInputMap2.CoordSystems=DataInputMap.CoordSystems;
DataInputMap2.xpts=new_xpts;
DataInputMap2.ypts=new_ypts;
DataInputMap2.max_pats=new_xpts*new_ypts;

DataInputMap2.PMap=reshape(DataInputMap.PMap(linears),new_ypts,new_xpts);
DataInputMap2.PMap=reshape(DataInputMap2.PMap',[],1);
%[DataInputMap2.YBeam_Map,DataInputMap2.XBeam_Map]=meshgrid([1:1:new_ypts],[1:1:new_xpts]);

%% Remap X-beam and Y-beam
original_x=DataInputMap2.XBeam_Map;
remapped_x=original_x;

original_y=DataInputMap2.YBeam_Map;
remapped_y=original_y;

x_mapping=[1:1:new_xpts];
y_mapping=[1:1:new_ypts];

for i=1:length(x_mapping)
    remapped_x(original_x==col_locs(i))=x_mapping(i);
end

for i=1:length(y_mapping)
    remapped_y(original_y==row_locs(i))=y_mapping(i);
end

DataInputMap2.XBeam=remapped_x;
DataInputMap2.YBeam=remapped_y;


%%
%DataInputMap2.PMap=DataInputMap2.PMap2;
% DataInputMap2.PMap=reshape([1:1:DataInputMap2.max_pats],new_ypts,new_xpts);
% DataInputMap2.PMap=reshape(DataInputMap2.PMap',1,[])';

%% MicroscopeData2
MicroscopeData2=MicroscopeData;
MicroscopeData2.NCOLS=new_xpts;
MicroscopeData2.NROWS=new_ypts;
MicroscopeData2.XSTEP=factor*MicroscopeData.XSTEP;
MicroscopeData2.YSTEP=factor*MicroscopeData.YSTEP;
MicroscopeData2.NPoints=DataInputMap2.max_pats;
MicroscopeData2.SEPixelSizeX=MicroscopeData.SEPixelSizeX*factor;
MicroscopeData2.SEPixelSizeY=MicroscopeData.SEPixelSizeY*factor;
MicroscopeData.SEMImage_T=MicroscopeData.SEMImage';
MicroscopeData2.SEMImage=reshape(MicroscopeData.SEMImage_T(linears)',new_ypts,new_xpts)';
%MicroscopeData2.SEMImage=reshape(MicroscopeData2.SEMImage,[],1);

%% write to hdf5
OutputUser=InputUser;

%new variable names
OutputUser.HDF5_file=[OutputUser.HDF5_file(1:end-3),'_x',num2str(factor),'_spatialreshape.h5'];
%we assume that this is an h5 file, hdf5 will cause issues with the (end-3)
OutputUser.DataName=InputUser.HDF5_file(1:end-3);
OutputUser.HDF5FullFile=fullfile(OutputUser.HDF5_folder,OutputUser.HDF5_file);

dtype='/EBSD/Data/'; %EBSD data location
htype='/EBSD/Header/'; %Header data location (e.g. microscope settings)

%% Map Data
%pattern centre
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,dtype,'DD',DataInputMap2.DD); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,dtype,'PCX',DataInputMap2.PCX); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,dtype,'PCY',DataInputMap2.PCY); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,dtype,'PMap',DataInputMap2.PMap);

%indexing data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,dtype,'MAD',DataInputMap2.MAD); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,dtype,'MADPhase',double(DataInputMap2.MADPhase)); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,dtype,'NIndexedBands',double(DataInputMap2.NIndexedBands)); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,dtype,'RadonQuality',DataInputMap2.RadonQuality); %write the data

%orientation data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,dtype,'PHI',DataInputMap2.PHI); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,dtype,'phi2',DataInputMap2.phi2); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,dtype,'phi1',DataInputMap2.phi1); %write the data

%beam data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,dtype,'X BEAM',double(DataInputMap2.XBeam-1)); %subtract 1 from the data because of the 0 indexing
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,dtype,'Y BEAM',double(DataInputMap2.YBeam-1)); %subtract 1 from the data because of the 0 indexing
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,dtype,'X SAMPLE',DataInputMap2.XSample); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,dtype,'Y SAMPLE',DataInputMap2.YSample); %write the data

%% Header Data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'MADMax',MicroscopeData2.MADMax); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'MaxRadonBandCount',MicroscopeData2.MaxRadonBandCount); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'CameraTilt',MicroscopeData2.CameraTilt); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'SampleTilt',MicroscopeData2.SampleTilt); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'TotalTilt',MicroscopeData2.TotalTilt); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'KV',MicroscopeData2.KV); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'NCOLS',MicroscopeData2.NCOLS); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'NROWS',MicroscopeData2.NROWS); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'XSTEP',MicroscopeData2.XSTEP); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'YSTEP',MicroscopeData2.YSTEP); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'NPoints',MicroscopeData2.NPoints); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'PatternHeight',MicroscopeData2.PatternHeight); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'PatternWidth',MicroscopeData2.PatternWidth); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'Magnification',MicroscopeData2.Magnification); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'MapStepFactor',MicroscopeData2.MapStepFactor); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'SEPixelSizeX',MicroscopeData2.SEPixelSizeX); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'SEPixelSizeY',MicroscopeData2.SEPixelSizeY); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'WD',MicroscopeData2.WD); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'SEM Image',MicroscopeData2.SEMImage); %write the data


%% Write data
dtype='uint16';

EBSPData.HDF5_loc=fullfile(OutputUser.HDF5_folder,OutputUser.HDF5_file);
EBSPData.numpats=MicroscopeData2.NPoints;

h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'PatternHeight',MicroscopeData2.PatternHeight);
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'PatternWidth',MicroscopeData2.PatternWidth);

%%
disp('Pattern conversion started');
for p=1:MicroscopeData2.NPoints
    
    pattern_number=DataInputMap2.PMap(p);
    new_pat=bReadEBSP(EBSPData_original,pattern_number);
    [EDSData_cor,EDSData_raw] = bReadEDX(EBSPData_original,pattern_number,channum);
    
    q=p;
    
    if p == 1
        EBSPData.made=0;
        bWriteEBSP(EBSPData,p,uint16(new_pat),dtype);
        bWriteSpectrum(EBSPData,p,EDSData_raw,EDSData_cor,channum,dtype);
    else
        EBSPData.made=1;
        bWriteEBSP(EBSPData,p,uint16(new_pat),dtype);
        bWriteSpectrum(EBSPData,p,EDSData_raw,EDSData_cor,channum,'double');
    end
end

%% Read the data back - warnings will tell you if things have failed
[ DataInputMap2_new,MicroscopeData2_new,PhaseData,EBSPData ] = bReadHDF5( OutputUser );
%check that we can read this
[ EBSDPat_back ] = bReadEBSP(EBSPData,1);
[ EDSSpec_back ] = bReadEDX(EBSPData,1,channum);

%%
InputUser.HDF5_folder='I:\TomMcA\CarbideSet'; %folder where the h5/bcf is located
InputUser.HDF5_file='LR12ht_8kX_1.h5';
[MapData,MicroscopeData,PhaseData,EBSPData ]=bReadHDF5( InputUser );
[DataInputMap] = EBSD_Map(MapData,MicroscopeData);

InputUser_comp.HDF5_folder='I:\TomMcA\CarbideSet'; %folder where the h5/bcf is located
InputUser_comp.HDF5_file='LR12ht_8kX_1_x3_spatialreshape.h5';
[MapData2,MicroscopeData2,PhaseData2,EBSPData2 ]=bReadHDF5( InputUser_comp );
[DataInputMap2] = EBSD_Map(MapData2,MicroscopeData2);

imagesc(MicroscopeData2.SEMImage)


% This code is copyright Alex Foden and Ben Britton 09/04/2019
% Do not distribute.
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
% Modified by TPM 20/02/2020

function [RTM]= RTM_run(InputUser,RTM,Settings_CorX,MicroscopeData,MapInfo,time1)

%%
screensize=RTM.screensize;
Sampling_Freq=RTM.Sampling_Freq;
iterations=RTM.iterations;
LPTsize=RTM.LPTsize;
Bin_loc=RTM.Bin_loc;
Phase_Folder=RTM.Phase_Folder;

PatternInfo.ScreenWidth=RTM.screensize;
PatternInfo.ScreenHeight=RTM.screensize;

if isfield(RTM,'InputPats')==1
    InputPats=RTM.InputPats;
end

if isfield(RTM,'ComponentMap')==1
    ComponentMap=RTM.ComponentMap;
end

if isfield(MapInfo,'Refined')==0
    MapInfo.Refined=0;
end

%Initialise some more stuff...
%This complication is necessary to provide generality to this function - it
%can either perform RTM on exisitng input patterns, or load them directly
%from a specified file.

if exist('InputPats')==1
    skipRead=1;
    exist_inputpats=1;
else
    skipRead=0;
    exist_inputpats=0;
    InputPats=0; %to let the parfor compile
end

if exist('ComponentMap')==1
    exist_components=1;
else
    exist_components=0;
    ComponentMap=0; %to let the parfor compile
end

%Set up the FFT filters for FFT and LPT
[ SettingsXCF, correction, SettingsXCF2 ] = FFT_Filter_settings( screensize, LPTsize );

%Define all rotation matrices needed in the code
RTM.Rz=@(theta)[cos(theta) sin(theta) 0;-sin(theta) cos(theta) 0;0 0 1]; %z rotation
RTM.Rx=@(theta)[1 0 0;0 cos(theta) sin(theta);0 -sin(theta) cos(theta)]; %x rotation
RTM.Ry=@(theta)[cos(theta) 0 sin(theta);0 1 0; -sin(theta) 0 cos(theta)]; %y rotation

InputUser.BCF_folder=InputUser.HDF5_folder; %this is used for the BCF reader

%update the pattern centre for square cropping
%If the PC has already been refined this isn't necessary.
if MapInfo.Refined==0
    pTime('Updating the PC map to adjust for square data',time1);
    [ Data_InputMap ] = PC_square( MapInfo.EBSPData, MapInfo.Data_InputMap,Settings_CorX );
else
    pTime('Matching with refined PC map',time1);
    Data_InputMap=MapInfo.Data_InputMap;
end
    
%Define the tilt of the detector
Detector_tilt = RTM.Rx(MicroscopeData.TotalTilt);

%% Handle whether to load patterns as we go, or use an input stack, and preallocate the FFT variables
if exist_inputpats==0
    P_list=Data_InputMap.PMap(:);
    tot_P_list=numel(P_list);
else
    tot_P_list=size(InputPats,3);
end
    
%Get the average PC - needed for geometry
Mean_PC=[mean(Data_InputMap.PCX(:),'omitnan'),mean(Data_InputMap.PCY(:),'omitnan'),mean(Data_InputMap.DD(:),'omitnan')];
[ Mean_EBSD_geom ] = EBSP_Gnom( PatternInfo,Mean_PC);

%need to do this to the size of an array
[~,XCF_data_fill] = refine_prep(zeros(Settings_CorX.size,Settings_CorX.size),SettingsXCF,RTM);
Full_PC=zeros(tot_P_list,3);
%preallocate the transform structure
Pat_Ref_r = struct('FFT', zeros(numel(XCF_data_fill), numel(XCF_data_fill)), 'logp', zeros(RTM.LPTsize,RTM.LPTsize));

%% Set up the matching...

Full_PC=zeros(tot_P_list,3);
PMap=MapInfo.Data_InputMap.PMap;

PCX=Data_InputMap.PCX;
PCY=Data_InputMap.PCY;
DD=Data_InputMap.DD;

if exist_inputpats==0
    
    %Load ALL the patterns in the map...
    parfor P_test=1:tot_P_list
        [P_r,P_c]=fPointFind(P_list(P_test),PMap);
        %read the EBSP
        [ Pat_Ref ] = bReadEBSP(MapInfo.EBSPData,P_list(P_test));
        %correct this pattern
        [ Pat_Ref_BG] = EBSP_BGCor( Pat_Ref,Settings_CorX );
        
        %generate the correct geometry for the full map, as we can
        Full_PC(P_test,:)=[PCX(P_r,P_c),PCY(P_r,P_c),DD(P_r,P_c)];
        
        %extract the LPT and the FFT
        [Pat_Ref_r(P_test)] = refine_prep(Pat_Ref_BG,SettingsXCF,RTM);
    end

else %We're loading from an input stack of patterns
    parfor P_test=1:tot_P_list
        %Take the EBSP from the input stack
        Pat_Ref_BG=InputPats(:,:,P_test);
        %Get the geometry
        
        %do we have a component map we might as well use?
       if exist_inputpats==1 && exist_components==1
            locs=find(ComponentMap==P_test);
            Full_PC(P_test,:) = [median(PCX(locs)),median(PCY(locs)),median(DD(locs))];

       else % get the median from the entire map
            Full_PC(P_test,:) = [median(PCX(:),'omitnan'),median(PCY(:),'omitnan'),median(DD(:),'omitnan')];
       end   
    
    %extract the LPT and the FFT
    [Pat_Ref_r(P_test)] = refine_prep(Pat_Ref_BG,SettingsXCF,RTM);
    end

end

%% Run the matching...
%pTime(['Testing against the libraries'],t1);

num_phases=length(InputUser.Phases);

%preallocate
G_Out=zeros(3,3,tot_P_list,num_phases);
PH_Out=zeros(tot_P_list,num_phases);
EAng_Out=zeros(3,tot_P_list,num_phases);

for num_lib=1:num_phases %run through each phase

    PhaseInput=InputUser.Phases{num_lib};
    pTime(['Testing against ' PhaseInput],time1);

    %build the phase
    [ ~,~,~,~,~, RTM_info ] = Phase_Builder_RTM( {PhaseInput},RTM.Phase_Folder);
    cs_phase=loadCIF(RTM_info.cif_file);

    %generate the simulated pattern cube
    [screen_int] = Cube_Generate(RTM_info.bin_file,RTM_info.isHex);

    %populate SO(3)
    [ library_G ] = SO3_rotmat_gen( cs_phase,RTM.Sampling_Freq);

    %generate the FFT of all patterns for the SO(3) library usin mean PC
    [ template_library ] = Library_Gen(Mean_EBSD_geom,screen_int,RTM_info.isHex,library_G,2,SettingsXCF);

    % get adjustment factor for radial crop
    radial_adjustment=fPHFac(Mean_EBSD_geom,Settings_CorX,screen_int,RTM,SettingsXCF);
    
    
    switch RTM.parsearch % library serach in parallel: RTM.parsearch=1
        case 1 %run through each pattern in series, with library search in parallel
            
            for P_test=1:tot_P_list
                %test the library
                [G_SO3,Library_PH]=fLibraryTest(template_library,library_G,Pat_Ref_r(P_test).FFT,SettingsXCF,XCF_data_fill,1);

                %generate the geometry for this location
                [ EBSD_geom ] = EBSP_Gnom( PatternInfo,Full_PC(P_test,:));
                %refine
                [G_Refined,regout_R] = refine5(Pat_Ref_r(P_test),EBSD_geom,Full_PC(P_test,:),G_SO3,SettingsXCF,screen_int,RTM_info.isHex,RTM);

                %store the orientation
                G_Out(:,:,P_test,num_lib)=G_Refined;
                %store the PH
                PH_Out(P_test,num_lib)=regout_R(4)./radial_adjustment;
            end
            
        case 2 % analyse PATTERNS in parallel, with lib search in series
            
            parfor P_test=1:tot_P_list
                %test the library
                [G_SO3,Library_PH]=fLibraryTest(template_library,library_G,Pat_Ref_r(P_test).FFT,SettingsXCF,XCF_data_fill,0);

                %generate the geometry for this location
                [ EBSD_geom ] = EBSP_Gnom( PatternInfo,Full_PC(P_test,:));
                %refine
                [G_Refined,regout_R] = refine5(Pat_Ref_r(P_test),EBSD_geom,Full_PC(P_test,:),G_SO3,SettingsXCF,screen_int,RTM_info.isHex,RTM);

                %store the orientation
                G_Out(:,:,P_test,num_lib)=G_Refined;
                %store the PH
                PH_Out(P_test,num_lib)=regout_R(4)./radial_adjustment;
            end
            
        case 0 %do everything in series
            
            for P_test=1:tot_P_list
                %test the library
                [G_SO3,Library_PH]=fLibraryTest(template_library,library_G,Pat_Ref_r(P_test).FFT,SettingsXCF,XCF_data_fill,0);

                %generate the geometry for this location
                [ EBSD_geom ] = EBSP_Gnom( PatternInfo,Full_PC(P_test,:));
                %refine
                [G_Refined,regout_R] = refine5(Pat_Ref_r(P_test),EBSD_geom,Full_PC(P_test,:),G_SO3,SettingsXCF,screen_int,RTM_info.isHex,RTM);

                %store the orientation
                G_Out(:,:,P_test,num_lib)=G_Refined;
                %store the PH
                PH_Out(P_test,num_lib)=regout_R(4)./radial_adjustment;
                
                
                
            end
            
    end %end the switch

    %convert G_Out to Euler angles - seperate in case we want to
    %parallelise above loop.
    for P_test=1:tot_P_list
        EAng_Out(:,P_test,num_lib)=conv_G_to_EA(G_Out(:,:,P_test,num_lib)*inv(Detector_tilt));
    end

end
    

% Or match HERE if input patterns are provided
%     
%     PCXs=Data_InputMap.PCX;
%     PCYs=Data_InputMap.PCY;
%     DDs=Data_InputMap.DD;
%     PMaps=Data_InputMap.PMap;
% 
%     PatternInfo.ScreenWidth=screensize;
%     PatternInfo.ScreenHeight=screensize;
%     
%     
%    parfor pat=1:testlib_size
% 
%        if exist_inputpats==1 && exist_components==1
% 
%             PatternCor=InputPats(:,:,pat);
%             locs=find(ComponentMap==pat);
%             PC_X = median(PCXs(locs));
%             PC_Y = median(PCYs(locs));
%             PC_Z = median(DDs(locs));
%             pnum=median(PMaps(locs));
%             PC_pat=[PC_X,PC_Y,PC_Z];
% 
%        else
% 
%             PatternCor=InputPats(:,:,pat);
%             PC_X = median(PCXs);
%             PC_Y = median(PCYs);
%             PC_Z = median(DDs);
%             pnum=median(PMaps);
%             PC_pat=[PC_X,PC_Y,PC_Z];
% 
%         end
%         
%         
%         %build the geometry
%         [ EBSP_pat ] = EBSP_Gnom( PatternInfo,PC_pat);
% 
%         %match from library
%         [rotmat_1st, shift_xFFT, shift_yFFT,rotmat_Peakheight] = index_FFT2(SettingsXCF, PatternCor, library_G, template_library_par.Value);
% 
%         %refine the orientation
%         [rotmat_best,regout] = refine4(PatternCor,PC_av,EBSP_pat,rotmat_1st,Rx,Ry,Rz,SettingsXCF,SettingsXCF2,correction,screen_int_par.Value,RTM_info.isHex,LPTsize,iterations); 
% 
%         %convert to euler angles
%         eang_template=conv_G_to_EA(rotmat_best*inv(Detector_tilt));
% 
%         TemData(:,pat)=[eang_template regout]';
%     end
%     


% RTM_peakheight(phasenum,:)=TemData(7,:);
% Euler1(phasenum,:)=TemData(1,:);
% Euler2(phasenum,:)=TemData(2,:);
% Euler3(phasenum,:)=TemData(3,:);

%%

% RTM.Output.PeakHeight=RTM_peakheight;
% RTM.Output.euler1=Euler1;
% RTM.Output.euler2=Euler2;
% RTM.Output.euler3=Euler3;

%%

RTM.Output.PeakHeight=PH_Out;
RTM.Output.Eulers=EAng_Out;

%%
pTime('finished all phases', time1);
end

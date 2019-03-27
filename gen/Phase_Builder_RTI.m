function [ Crystal_UCell,Crystal_Family,Crystal_LUT,Settings_LUT,num_Phases, RTI_info ] = Phase_Builder_RTI( Phases,Phase_Folder, Bin_loc )
%PHASE_BUILD Build all the phases needed for the EBSD data
%Phases = a structure of phase file names
%Phase_Folder = folder where they exist

%% Versioning
%v1 - TBB 14/04/2017

%% run code
num_Phases=size(Phases,2);

for num_P=1:num_Phases
    %read the phases
    Phase_File=fullfile(Phase_Folder,[Phases{num_P} '.pha']);
    if ~exist(Phase_File) %#ok<EXIST>
        error(['This phase file does not exist: ' Phase_File])
    end
    
    %generate the unit cells
    [Crystal_UCell{num_P}]=Build_UCell(Phase_File); %#ok<*AGROW>
    %generate the plane families
    [Crystal_Family{num_P} ] = Build_Reflec( Crystal_UCell{num_P} );
    %generate the look up tables (LUTs)
    [Crystal_LUT{num_P}] = Build_LUT(Crystal_Family{num_P});
    %threshold for the LUT generation
    Settings_LUT{num_P}.thresh_trig=1/Crystal_UCell{num_P}.efac;
end

%% RTI extra information needed
delimiter = '';
formatSpec = '%s%[^\n\r]';
fileID = fopen(Phase_File,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);
fclose(fileID);
PhaseRaw = [dataArray{1:end-1}];

array_num=1:size(PhaseRaw,1);

%Extract the cif file
ix=array_num(contains(PhaseRaw,'$cif','IgnoreCase',true));
RTI_info.cif_file=PhaseRaw{ix+1};

%Extract the bin file name
ix=array_num(contains(PhaseRaw,'$dynamic','IgnoreCase',true));
bin_file_name=PhaseRaw{ix+1};
RTI_info.bin_file = [Bin_loc '\' bin_file_name '.bin'];

%Extract the ishex information
ix=array_num(contains(PhaseRaw,'$isHex','IgnoreCase',true));
RTI_info.isHex=PhaseRaw{ix+1};

end


function [ Crystal_UCell,Crystal_Family,Crystal_LUT,Settings_LUT,num_Phases ] = Phase_Builder( Phases,Phase_Folder )
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

end


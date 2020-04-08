function [screen_int,RTM_info,cslist] = Build_Interps(InputUser)
%BUILD_INTERPS Creates interpolants and phase lists for phases
%
% [screen_int,RTM_info,cslist] = Build_Interps(InputUser,RTM)
% 
% INPUTS
% InputUser.Phases = {'Ferrite','Austenite'} 
%                     - All phases as needed
% InputUser.Phase_Folder = 'E:\Ben\phases\' 
%                           - location where the phase folder is found
%
% OUTPUTS
% screen_int{n} = screen_int interpolants
% RTM_info{n} = structure with cif file, bin file, and isHex
% cslist = cslist for these phases
%
% Created TBB 07/04/2020

cslist={'notIndexed'};
if ~isfield(InputUser,'Phases')
    error('InputUser.Phases not created')
end

num_phases=numel(InputUser.Phases);
for n=1:num_phases
    PhaseInput=InputUser.Phases{n};
    [ ~,~,~,~,~, RTM_info(n) ] = Phase_Builder_RTM( {PhaseInput},InputUser.Phase_Folder);
    cs_phase=loadCIF(RTM_info(n).cif_file);
    cslist{end+1}=cs_phase;
    [screen_int(n)] = Cube_Generate(RTM_info(n).bin_file,RTM_info(n).isHex);
end


end


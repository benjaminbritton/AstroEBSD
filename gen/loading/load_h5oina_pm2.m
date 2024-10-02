function [ebsd_original,dataset_header,ebsd_patternmatched,h5_original,h5_patternmatch]=load_h5oina_pm2(fname_slice_data,fname_folder,warningOn)
% LOAD_H5OINA_PM2 Loads h5oina data, and pattern matched data
% Option 1 - pattern matched data should be in the same folder and with the file name of the original
%            but an additional '-PatternMatching' appended
% Option 2 - pattern matcched data is saved within the same file
%
% INPUTS
%   fname_slice_data    - h5oina file name  (excluding the folder and h5oina extension)
%   fname_folder        - folder location of the h5oina file
%   WarningOn           - output warnings
%
% OUTPUTS
%   ebsd_original       - ebsd container (original data)
%   dataset_header      - header
%   ebsd_patternmatched - ebsd container (pattern matched)
%   h5_original         - file name original file path
%   h5_patternmatch     - file name pattern matched file path
%
% Authors: TBB & RB
% Last update: 03/09/2024 RB
%% Setup

%do not change unless you really need to
fset_match_suffix='-PatternMatching';
fset_fileformat='.h5oina';

%h5 full files
h5_original=fullfile(fname_folder,[fname_slice_data fset_fileformat]);
h5_patternmatch=fullfile(fname_folder,[fname_slice_data fset_match_suffix fset_fileformat]);

%check file exists
if exist(h5_original,'file') ~= 2
    if warningOn==1
        error([h5_original ' can not be found'])
    end
end

%check if extra pattern matched file exists
if exist(h5_patternmatch,'file') == 2
    pm_data=1;
else
    if warningOn==1
        warning([h5_patternmatch ' can not be found - ignore if data not pattern matched'])
    end
    pm_data=0;
end

%% Load data

% load data - main file (inc. pattern matched if present)
[ebsd_o,EBSDdata,opt_o,ebsdpm2,EBSDdata2]= loadEBSD_h5oina_getopt_i(h5_original,warningOn);

% load pattern matched data if in separate file
if ~isempty(fieldnames(ebsdpm2))
    ebsd_pm=ebsdpm2;
elseif pm_data == 1 % pattern matched (2 file version)
    [ebsd_pm,~,~,~] = loadEBSD_h5oina_PatMatch_i(h5_patternmatch,ebsd_o,opt_o,warningOn) ;
end

%% create header

% header data
dataset_header=ebsd_o.opt.Header;

% oiginal ebsd
ebsd_original=ebsd_o;

% pattern matched ebsd
if pm_data == 1 % 2 file pattern matched
    ebsd_patternmatched=ebsd_pm;
    disp([fname_slice_data ' - ebsd data loaded including pattern matched data'])
elseif ~isempty(fieldnames(ebsdpm2)) % 1 file pattern matched
    ebsd_patternmatched=ebsdpm2;
    disp([fname_slice_data ' - ebsd data loaded including pattern matched data'])
else % no pattern matching
    ebsd_patternmatched=ebsd_original;
    disp([fname_slice_data ' - ebsd data loaded'])
end

end

function [ebsd,EBSDdata,opt,ebsdpm,EBSDdata2] = loadEBSD_h5oina_getopt_i(fname,warningOn,varargin)
% read HKL *.h5oina hdf5 file
% documented here: https://github.com/oinanoanalysis/h5oina/blob/master/H5OINAFile.md
% note that Matlab < R2021b does not handle hdf5 v1.10 and one needs to use hdf5format_convert
% (https://github.com/HDFGroup/hdf5) on the input file to prevent Matlab from fatally crashing

% TODO
% 1) Test if EBSDheader.Specimen_Orientation_Euler does what it's supposed
%    to do -> see below
% 2) find a solution if multiple ebsd datasets are contained, export to a
%    cell?
% 3) decide what header data to use and how to display it? Fix display for
% the header to be shown correctly (bc. ebsd.opt.Header sort of works)


%create variables in case they are needed
ebsdpm=struct;
EBSDdata2=struct;
opt2=struct;


all = h5info(fname);

% task: find all groups called EBSD, therein header and PHASE
EBSD_index = {};
EDS_index = {};
PM_index = {};

n = 1;
m = 1;
p = 1;

%search for EBSD/EDS/PM data
for i = 1:length(all.Groups) % map site on sample
    if ~isempty(all.Groups(i).Groups) % data on map site ('EBSD, EDS, Electron iamge etc)
        for j=1:length(all.Groups(i).Groups)
            if contains(all.Groups(i).Groups(j).Name,'EBSD')
                EBSD_index{n} = [i j];
                n = n+1;
            end
            if contains(all.Groups(i).Groups(j).Name,'EDS')
                EDS_index{m} = [i j];
                m = m+1;
            end

            if contains(all.Groups(i).Groups(j).Name,'Data Processing')
                PM_index{p} = [i j];
                p = p+1;
            end

        end

    end
end

%EBSD data reading
if length(EBSD_index) > 1
    disp('more than 1 EBSD dataset in the file, output will be a cell')
end

for k = 1 :length(EBSD_index) % TODO: find a good way to write out multiple datasets
    %EBSD dataset
    EBSD_data = all.Groups(EBSD_index{k}(1)).Groups(EBSD_index{k}(2)).Groups(1);

    %EBSD header
    EBSD_header = all.Groups(EBSD_index{k}(1)).Groups(EBSD_index{k}(2)).Groups(2);


    if ~isempty(EDS_index) & EDS_index{k}(1) == EBSD_index{k}(1)

        % EDS times and coordiantes - not used for now
        % eds_tc = all.Groups(EDS_index{k}(1)).Groups(EDS_index{k}(2)).Groups(1)

        % EDS header
        EDS_header = all.Groups(EDS_index{k}(1)).Groups(EDS_index{k}(2)).Groups(2);

        % EDS data
        EDS_data= all.Groups(EDS_index{k}(1)).Groups(EDS_index{k}(2)).Groups(1).Groups;

    end


    % read all EBSD data
    EBSDdata = struct;
    for thing = 1:length(EBSD_data.Datasets)
        sane_name = regexprep(EBSD_data.Datasets(thing).Name,' |-|,|:|%|~|#','_');
        if strcmpi(sane_name,'Processed_Patterns')
            if warningOn ==1
                disp('patterns present but not loaded')
            end
        else
            EBSDdata.(sane_name)=double(h5read(fname,[EBSD_data.Name '/' EBSD_data.Datasets(thing).Name]));
        end
    end

    %read EBSD header
    EBSDheader = struct;
    for thing = 1:length(EBSD_header.Datasets)
        sane_name = regexprep(EBSD_header.Datasets(thing).Name,' |-|,|:|%|~|#','_');
        content = h5read(fname,[EBSD_header.Name '/' EBSD_header.Datasets(thing).Name]);
        if any(size(content) ~=1) & isnumeric(content)
            content = reshape(content,1,[]);
        end
        EBSDheader.(sane_name) = content;
    end



    %read the EDS data
    if ~isempty(EDS_index) & EDS_index{k}(1) == EBSD_index{k}(1)
        try
            %read EDS data
            EDSdata = struct;
            for thing = 1:length(EDS_data.Datasets)
                sane_name = regexprep(EDS_data.Datasets(thing).Name,' |-|,|:|%|~|#','_');
                EDSdata.(sane_name)=double(h5read(fname,[EDS_data.Name '/' EDS_data.Datasets(thing).Name]));
            end
            %read EDS header
            EDSheader = struct;
            for thing = 1:length(EDS_header.Datasets)
                sane_name = regexprep(EDS_header.Datasets(thing).Name,' |-|,|:|%|~|#','_');
                content = h5read(fname,[EDS_header.Name '/' EDS_header.Datasets(thing).Name]);
                if any(size(content) ~=1) & isnumeric(content)
                    content = reshape(content,1,[]);
                end
                EDSheader.(sane_name) = content;
            end
        catch
            if warningOn==1
                warning('Error loading EDS data - catch loop invoked')
            end
        end
    end

    EBSDphases = struct;
    phases = all.Groups(EBSD_index{k}(1)).Groups(EBSD_index{k}(2)).Groups(2).Groups(1);
    %   ----------------

    CS{1}='notIndexed';
    for phaseN = 1:length(phases.Groups)

        phase_name=phases.Groups(phaseN).Name;
        fslash_loc=strfind(phase_name,'/');
        phase_num=str2num(phase_name(fslash_loc(end)+1:end));

        pN = ['phase_' num2str(phaseN)];




        EBSDphases.(pN)= struct;
        for j = 1:length(phases.Groups(phaseN).Datasets)
            sane_name = regexprep(phases.Groups(phaseN).Datasets(j).Name,' |-|,|:|%|~|#','_');
            content = h5read(fname,[phases.Groups(phaseN).Name '/' phases.Groups(phaseN).Datasets(j).Name]);
            EBSDphases.(pN).(sane_name) = content;
        end

        %remove underscores
        phase_name_no_und=char(EBSDphases.(pN).Phase_Name);
        score_loc=strfind(phase_name_no_und,'_');
        phase_name_no_und(score_loc)='-';

        % the angle comes in single precision. make sure something
        % sufficiently close to 90 resp. 120 does not end up with
        % rounding errors instead of using the 'force' option

        langle = double(EBSDphases.(pN).Lattice_Angles');
        % if pN=='phase_1'
        %     EBSDphases.(pN).Space_Group=192;
        % else
        % end
        csm = crystalSymmetry('SpaceId',EBSDphases.(pN).Space_Group);
        if strcmp(csm.lattice,'trigonal') | strcmp(csm.lattice,'hexagonal')
            langle(isnull(langle-2/3*pi,1e-7))=2/3*pi;
        else
            langle(isnull(langle-pi/2,1e-7))=pi/2;
        end

        CS{phase_num} = crystalSymmetry('SpaceId',EBSDphases.(pN).Space_Group, ...
            double(EBSDphases.(pN).Lattice_Dimensions'),...
            'Mineral',phase_name_no_und);%%             langle,...
        %             'X||a*','Y||b', 'Z||C');
    end


    % write out first EBSD dataset
    % EBSDheader.Specimen_Orientation_Euler: this should be the convention to map
    % CS1 (sample surface) into CS0 (sample primary),
    % CS2 into CS1 should be absolute orientation
    % TODO! make sure those rotations are correctly applied, possibly
    % EBSDheader.Specimen_Orientation_Euler

    rc = rotation.byEuler(double(EBSDheader.Specimen_Orientation_Euler*degree)); % what definition? Bunge?

    % set up EBSD data
    rot = rc*rotation.byEuler(EBSDdata.Euler');
    phase = EBSDdata.Phase;
    opt=struct;

    % read some fields
    EBSD_fieldnames=fieldnames(EBSDdata);

    num_fields=size(EBSD_fieldnames,1);
    for n = 1: num_fields
        s=EBSDdata.(EBSD_fieldnames{n});
        if size(s,2) == 1 && size(s,1) == numel(phase)
            try
                opt.(EBSD_fieldnames{n})=s;
            catch
            end

        end
    end

    %not put the Euler angles in the options too
    if isfield(EBSDdata,'Euler')
        opt.euler1=EBSDdata.Euler(1,:);
        opt.euler2=EBSDdata.Euler(2,:);
        opt.euler3=EBSDdata.Euler(3,:);
    end

    try
        opt.x=opt.Beam_Position_X;
        opt.y=opt.Beam_Position_Y;
    catch
        error('Beam Positions Not Loaded - was this a pattern match file?')
    end


    try
        % read all PM data
        PMdata = struct;
        for thing = 1:length(EBSD_data.Datasets)
            sane_name = regexprep(EBSD_data.Datasets(thing).Name,' |-|,|:|%|~|#','_');
            if strcmpi(sane_name,'Processed_Patterns')
                if warningOn==1
                    disp('patterns present but not loaded')
                end
            else
                PMdata.(sane_name)=double(h5read(fname,[EBSD_data.Name '/' EBSD_data.Datasets(thing).Name]));
            end
        end

        %read the pattern matching data if its there
        for k = 1 :length(PM_index) % TODO: find a good way to write out multiple datasets
            num_groups=size(all.Groups(PM_index{k}(1)).Groups(PM_index{k}(2)).Groups,1);
            for n=1:num_groups
                groupname=all.Groups(PM_index{k}(1)).Groups(PM_index{k}(1)).Groups(n).Name;
                fslashloc=strfind(groupname,'/');
                subgroup_end=groupname(fslashloc(end)+1:end);

                if strcmpi(subgroup_end,'Data')
                    PM_data=all.Groups(PM_index{k}(1)).Groups(PM_index{k}(1)).Groups(n);
                end

                if strcmpi(subgroup_end,'Header')
                    PM_header=all.Groups(PM_index{k}(1)).Groups(PM_index{k}(1)).Groups(n);
                end

                if strcmpi(subgroup_end,'Pattern Matching')
                    PM_PatMatch=all.Groups(PM_index{k}(1)).Groups(PM_index{k}(1)).Groups(n);

                    num_groups2=size(PM_PatMatch.Groups,1);

                    for n2=1:num_groups2
                        groupname=PM_PatMatch.Groups(n2).Name;
                        fslashloc=strfind(groupname,'/');
                        subgroup_end=groupname(fslashloc(end)+1:end);

                        if strcmpi(subgroup_end,'Data')
                            PM_PatMatch2=PM_PatMatch.Groups(n2);
                        end
                    end
                end
            end

        end

        %read PM header
        PMheader = struct;
        for thing = 1:length(PM_header.Datasets)
            sane_name = regexprep(PM_header.Datasets(thing).Name,' |-|,|:|%|~|#','_');
            content = h5read(fname,[PM_header.Name '/' PM_header.Datasets(thing).Name]);
            if any(size(content) ~=1) & isnumeric(content)
                content = reshape(content,1,[]);
            end
            PMheader.(sane_name) = content;
        end

        %read PM data
        PMdata = struct;
        for thing = 1:length(PM_data.Datasets)
            sane_name = regexprep(PM_data.Datasets(thing).Name,' |-|,|:|%|~|#','_');
            % content = h5read(fname,[PM_data.Name '/' PM_data.Datasets(thing).Name]);
            % if any(size(content) ~=1) & isnumeric(content)
            %     content = reshape(content,1,[]);
            % end
            PMdata.(sane_name)=double(h5read(fname,[PM_data.Name '/' PM_data.Datasets(thing).Name]));

            % PMdata.(sane_name) = content;
        end

        %read PM_PatMatch

        for thing = 1:length(PM_PatMatch2.Datasets)
            sane_name = regexprep(PM_PatMatch2.Datasets(thing).Name,' |-|,|:|%|~|#','_');
            % content = h5read(fname,[PM_PatMatch2.Name '/' PM_PatMatch2.Datasets(thing).Name]);
            % if any(size(content) ~=1) & isnumeric(content)
            %     content = reshape(content,1,[]);
            % end
            PMdata.(sane_name)=double(h5read(fname,[PM_PatMatch2.Name '/' PM_PatMatch2.Datasets(thing).Name]));
            % PMdata.(sane_name) = content;
        end

        %load the pattern matched phases
        PMphases = struct;
        phases = PM_header.Groups(1);
        %   ----------------

        CS_pm{1}='notIndexed';
        %sort the phase names into a proper list


        for phaseN = 1:length(phases.Groups)
            %deal with a phase list that is longer than 9
            phase_name=phases.Groups(phaseN).Name;
            fslash_loc=strfind(phase_name,'/');
            phase_num=str2num(phase_name(fslash_loc(end)+1:end));

            pN = ['phase_' num2str(phase_num)];




            PMphases.(pN)= struct;
            for j = 1:length(phases.Groups(phaseN).Datasets)
                sane_name = regexprep(phases.Groups(phaseN).Datasets(j).Name,' |-|,|:|%|~|#','_');
                content = h5read(fname,[phases.Groups(phaseN).Name '/' phases.Groups(phaseN).Datasets(j).Name]);
                PMphases.(pN).(sane_name) = content;
            end

            %remove underscores
            phase_name_no_und=char(PMphases.(pN).Phase_Name);
            score_loc=strfind(phase_name_no_und,'_');
            phase_name_no_und(score_loc)='-';

            % the angle comes in single precision. make sure something
            % sufficiently close to 90 resp. 120 does not end up with
            % rounding errors instead of using the 'force' option

            langle = double(PMphases.(pN).Lattice_Angles');
            try
                csm = crystalSymmetry('SpaceId',PMphases.(pN).Space_Group);
            catch
                csm = crystalSymmetry('SpaceId',1);
                PMphases.(pN).Space_Group=1;
                % if warningOn==1
                warning(['Space group for PM phase ' pN ' not loaded, and set to 1 - IPF colours will be wrong if not corrected']);
                % end
            end

            if strcmp(csm.lattice,'trigonal') | strcmp(csm.lattice,'hexagonal')
                langle(isnull(langle-2/3*pi,1e-7))=2/3*pi;
            else
                langle(isnull(langle-pi/2,1e-7))=pi/2;
            end


            CS_pm{phase_num} = crystalSymmetry('SpaceId',PMphases.(pN).Space_Group, ...
                double(PMphases.(pN).Lattice_Dimensions'),...
                'Mineral',phase_name_no_und);%%             langle,...
            %             'X||a*','Y||b', 'Z||C');
        end


        % write out first EBSD dataset
        % EBSDheader.Specimen_Orientation_Euler: this should be the convention to map
        % CS1 (sample surface) into CS0 (sample primary),
        % CS2 into CS1 should be absolute orientation
        % TODO! make sure those rotations are correctly applied, possibly
        % EBSDheader.Specimen_Orientation_Euler

        rc2 = rotation.byEuler(double(EBSDheader.Specimen_Orientation_Euler*degree)); % what definition? Bunge?

        % set up EBSD data
        rot2 = rc2*rotation.byEuler(PMdata.Euler');
        phase2 = PMdata.Phase;
        opt2=struct;

        % read some fields
        EBSD_fieldnames=fieldnames(EBSDdata);

        num_fields=size(EBSD_fieldnames,1);
        for n = 1: num_fields
            s=EBSDdata.(EBSD_fieldnames{n});
            if size(s,2) == 1 && size(s,1) == numel(phase)
                try
                    opt2.(EBSD_fieldnames{n})=s;
                catch
                end

            end
        end

        % read some fields
        PM_fieldnames=fieldnames(PMdata);

        num_fields=size(PM_fieldnames,1);
        for n = 1: num_fields
            s=PMdata.(PM_fieldnames{n});
            if size(s,2) == 1 && size(s,1) == numel(phase)
                try
                    opt2.(PM_fieldnames{n})=s;
                catch
                end

            end
        end



        %not put the Euler angles in the options too
        if isfield(EBSDdata,'Euler')
            opt.euler1=EBSDdata.Euler(1,:);
            opt.euler2=EBSDdata.Euler(2,:);
            opt.euler3=EBSDdata.Euler(3,:);

            opt2.euler1=EBSDdata.Euler(1,:);
            opt2.euler2=EBSDdata.Euler(2,:);
            opt2.euler3=EBSDdata.Euler(3,:);
        end

        if isfield(PMdata,'Euler')
            opt.euler1_P=PMdata.Euler(1,:);
            opt.euler2_P=PMdata.Euler(2,:);
            opt.euler3_P=PMdata.Euler(3,:);

            opt2.euler1_P=PMdata.Euler(1,:);
            opt2.euler2_P=PMdata.Euler(2,:);
            opt2.euler3_P=PMdata.Euler(3,:);
        end
    catch
        if warningOn==1
            warning('Did not load pattern matched data')
        end
    end

    % if available, add EDS data
    if exist('EDSdata','var')
        eds_names = fieldnames(EDSdata);
        for j =1 :length(eds_names)
            opt.(eds_names{j}) = EDSdata.(eds_names{j});
        end
    end

    ebsdtemp = EBSD(rot,phase,CS,opt,'unitCell',calcUnitCell([opt.x,opt.y]));
    ebsdtemp.opt.Header = EBSDheader;

    if length(EBSD_index) > 1
        ebsd{k} = ebsdtemp;
    else
        ebsd = ebsdtemp;
    end

    try %pattern matched data
        opt2.Phase=phase2;
        opt2.x=opt.x;
        opt2.y=opt.y;
        ebsdtemp2 = EBSD(rot2,opt2.Phase,CS_pm,opt2,'unitCell',calcUnitCell([opt.x,opt.y]));
        ebsdtemp2.opt.Header = EBSDheader;
        ebsdtemp2.opt.Header_PM = PMheader;

        if length(EBSD_index) > 1
            ebsdpm{k} = ebsdtemp2;
        else
            ebsdpm = ebsdtemp2;
        end
    catch
        if warningOn==1
            warning('Pattern matched data not loaded')
        end
    end



end

end

function [ebsdDP,EBSDdata,DPdata1,DPdata3] = loadEBSD_h5oina_PatMatch_i(fname,ebsdOrig,optOrig,warningOn, varargin)
% read HKL *.h5oina hdf5 file
% documented here: https://github.com/oinanoanalysis/h5oina/blob/master/H5OINAFile.md
% note that Matlab < R2021b does not handle hdf5 v1.10 and one needs to use hdf5format_convert
% (https://github.com/HDFGroup/hdf5) on the input file to prevent Matlab from fatally crashing

% TODO
% 1) Test if EBSDheader.Specimen_Orientation_Euler does what it's supposed
%    to do -> see below
% 2) find a solution if multiple ebsd datasets are contained, export to a
%    cell?
% 3) decide what header data to use and how to display it? Fix display for
% the header to be shown correctly (bc. ebsd.opt.Header sort of works)

all = h5info(fname);

% task: find all groups called EBSD, therein header and PHASE
EBSD_index = {};
EDS_index = {};
DP_index = {}; % data processing parts

n = 1;
m=1;
%search for EBSD data
for i = 1:length(all.Groups) % map site on sample
    if ~isempty(all.Groups(i).Groups) % data on map site ('EBSD, EDS, Electron iamge etc)
        for j=1:length(all.Groups(i).Groups)
            if contains(all.Groups(i).Groups(j).Name,'EBSD')
                EBSD_index{n} = [i j];
                n = n+1;
            end
            if contains(all.Groups(i).Groups(j).Name,'EDS')
                EDS_index{m} = [i j];
                m = m+1;
            end
            if contains(all.Groups(i).Groups(j).Name,'Data Processing')
                DP_index{m} = [i j];
                m = m+1;
            end

        end

    end
end

if length(EBSD_index) > 1
    disp('more than 1 EBSD dataset in the file, output will be a cell')
end
for k = 1 :length(EBSD_index) % TODO: find a good way to write out multiple datasets
    %EBSD dataset
    EBSD_data = all.Groups(EBSD_index{k}(1)).Groups(EBSD_index{k}(2)).Groups(1);

    %EBSD header
    EBSD_header = all.Groups(EBSD_index{k}(1)).Groups(EBSD_index{k}(2)).Groups(2);

    %Mapsweeper data (like EBSD)
    DP_data1 = all.Groups(DP_index{k}(2)).Groups(DP_index{k}(2)).Groups(2);

    %Mapsweeper data (like EBSD Header)
    DP_data2 = all.Groups(DP_index{k}(2)).Groups(DP_index{k}(2)).Groups(3);

    %Mapsweeper data (cross corr )
    try
        DP_data3 = all.Groups(DP_index{k}(2)).Groups(DP_index{k}(2)).Groups(4).Groups(1);
    catch
        DP_data3= all.Groups(DP_index{k}(2)).Groups(DP_index{k}(2)).Groups(3).Groups(3);
    end

    if ~isempty(EDS_index) & EDS_index{k}(1) == EBSD_index{k}(1)

        % EDS times and coordiantes - not used for now
        % eds_tc = all.Groups(EDS_index{k}(1)).Groups(EDS_index{k}(2)).Groups(1)

        % EDS header
        EDS_header = all.Groups(EDS_index{k}(1)).Groups(EDS_index{k}(2)).Groups(2);

        % EDS data
        EDS_data= all.Groups(EDS_index{k}(1)).Groups(EDS_index{k}(2)).Groups(1).Groups;

    end


    % read all EBSD data
    EBSDdata = struct;
    for thing = 1:length(EBSD_data.Datasets)
        sane_name = regexprep(EBSD_data.Datasets(thing).Name,' |-|,|:|%|~|#','_');
        if strcmpi(sane_name,'Processed_Patterns')
            if warningOn==1
                disp('patterns present but not loaded')
            end
        else
            EBSDdata.(sane_name)=double(h5read(fname,[EBSD_data.Name '/' EBSD_data.Datasets(thing).Name]));
        end
    end

    %read EBSD header
    EBSDheader = struct;
    for thing = 1:length(EBSD_header.Datasets)
        sane_name = regexprep(EBSD_header.Datasets(thing).Name,' |-|,|:|%|~|#','_');
        content = h5read(fname,[EBSD_header.Name '/' EBSD_header.Datasets(thing).Name]);
        if any(size(content) ~=1) & isnumeric(content)
            content = reshape(content,1,[]);
        end
        EBSDheader.(sane_name) = content;
    end

    % read all EBSD data from MS
    DPdata1 = struct;
    for thing = 1:length(DP_data1.Datasets)
        sane_name = regexprep(DP_data1.Datasets(thing).Name,' |-|,|:|%|~|#','_');
        if strcmpi(sane_name,'Processed_Patterns')
            disp('patterns present but not loaded')
        else
            DPdata1.(sane_name)=double(h5read(fname,[DP_data1.Name '/' DP_data1.Datasets(thing).Name]));
        end
    end

    % read all EBSD header data from MS
    DPdata2 = struct;
    for thing = 1:length(DP_data2.Datasets)
        sane_name = regexprep(DP_data2.Datasets(thing).Name,' |-|,|:|%|~|#','_');
        content = h5read(fname,[DP_data2.Name '/' DP_data2.Datasets(thing).Name]);
        if any(size(content) ~=1) & isnumeric(content)
            content = reshape(content,1,[]);
        end
        DPdata2.(sane_name) = content;
    end

    % read all EBSD header data from MS
    DPdata3 = struct;
    for thing = 1:length(DP_data3.Datasets)
        sane_name = regexprep(DP_data3.Datasets(thing).Name,' |-|,|:|%|~|#','_');
        content = h5read(fname,[DP_data3.Name '/' DP_data3.Datasets(thing).Name]);
        if any(size(content) ~=1) & isnumeric(content)
            content = reshape(content,1,[]);
        end
        DPdata3.(sane_name) = content;
    end

    if ~isempty(EDS_index) & EDS_index{k}(1) == EBSD_index{k}(1)
        %read EDS data
        EDSdata = struct;
        for thing = 1:length(EDS_data.Datasets)
            sane_name = regexprep(EDS_data.Datasets(thing).Name,' |-|,|:|%|~|#','_');
            EDSdata.(sane_name)=double(h5read(fname,[EDS_data.Name '/' EDS_data.Datasets(thing).Name]));
        end
        %read EDS header
        EDSheader = struct;
        for thing = 1:length(EDS_header.Datasets)
            sane_name = regexprep(EDS_header.Datasets(thing).Name,' |-|,|:|%|~|#','_');
            content = h5read(fname,[EDS_header.Name '/' EDS_header.Datasets(thing).Name]);
            if any(size(content) ~=1) & isnumeric(content)
                content = reshape(content,1,[]);
            end
            EDSheader.(sane_name) = content;
        end

    end

    EBSDphases = struct;
    phases = all.Groups(EBSD_index{k}(1)).Groups(EBSD_index{k}(2)).Groups(2).Groups(1);

    %   ----------------

    CS{1}='notIndexed';
    for phaseN = 1:length(phases.Groups)
        pN = ['phase_' num2str(phaseN)];
        EBSDphases.(pN)= struct;
        for j = 1:length(phases.Groups(phaseN).Datasets)
            sane_name = regexprep(phases.Groups(phaseN).Datasets(j).Name,' |-|,|:|%|~|#','_');
            content = h5read(fname,[phases.Groups(phaseN).Name '/' phases.Groups(phaseN).Datasets(j).Name]);
            EBSDphases.(pN).(sane_name) = content;
        end

        % the angle comes in single precision. make sure something
        % sufficiently close to 90 resp. 120 does not end up with
        % rounding errors instead of using the 'force' option

        langle = double(EBSDphases.(pN).Lattice_Angles');
        if pN=='phase_1'
            EBSDphases.(pN).Space_Group=192;
        else
        end
        csm = crystalSymmetry('SpaceId',EBSDphases.(pN).Space_Group);
        if strcmp(csm.lattice,'trigonal') | strcmp(csm.lattice,'hexagonal')
            langle(isnull(langle-2/3*pi,1e-7))=2/3*pi;
        else
            langle(isnull(langle-pi/2,1e-7))=pi/2;
        end

        CS{phaseN} = crystalSymmetry('SpaceId',EBSDphases.(pN).Space_Group, ...
            double(EBSDphases.(pN).Lattice_Dimensions'),...
            'Mineral',char(EBSDphases.(pN).Phase_Name));%%             langle,...
        %             'X||a*','Y||b', 'Z||C');
    end
    %%
    DPphases = struct;
    phases1 = all.Groups(DP_index{k}(1)).Groups(DP_index{k}(2)).Groups(3).Groups(1);

    CS{1}='notIndexed';
    for phaseN = 1:length(phases1.Groups)
        pN = ['phase_' num2str(phaseN)];
        DPphases.(pN)= struct;
        for j = 1:length(phases1.Groups(phaseN).Datasets)
            sane_name = regexprep(phases1.Groups(phaseN).Datasets(j).Name,' |-|,|:|%|~|#','_');
            content = h5read(fname,[phases1.Groups(phaseN).Name '/' phases1.Groups(phaseN).Datasets(j).Name]);
            DPphases.(pN).(sane_name) = content;
        end

        % the angle comes in single precision. make sure something
        % sufficiently close to 90 resp. 120 does not end up with
        % rounding errors instead of using the 'force' option

        langle = double(DPphases.(pN).Lattice_Angles');
        if pN=='phase_1'
            DPphases.(pN).Space_Group=192;
        else
        end
        csm = crystalSymmetry('SpaceId',DPphases.(pN).Space_Group);
        if strcmp(csm.lattice,'trigonal') | strcmp(csm.lattice,'hexagonal')
            langle(isnull(langle-2/3*pi,1e-7))=2/3*pi;
        else
            langle(isnull(langle-pi/2,1e-7))=pi/2;
        end

        CS{phaseN} = crystalSymmetry('SpaceId',DPphases.(pN).Space_Group, ...
            double(DPphases.(pN).Lattice_Dimensions'),...
            'Mineral',char(DPphases.(pN).Phase_Name));%%             langle,...
        %             'X||a*','Y||b', 'Z||C');
    end

    %%
    % write out first EBSD dataset
    % EBSDheader.Specimen_Orientation_Euler: this should be the convention to map
    % CS1 (sample surface) into CS0 (sample primary),
    % CS2 into CS1 should be absolute orientation
    % TODO! make sure those rotations are correctly applied, possibly
    % EBSDheader.Specimen_Orientation_Euler

    rc = rotation.byEuler(double(EBSDheader.Specimen_Orientation_Euler*degree)); % what definition? Bunge?

    % set up EBSD data
    rot = rc*rotation.byEuler(EBSDdata.Euler');
    phase = EBSDdata.Phase;
    opt=struct;

    % set up EBSD data
    rotDP = rc*rotation.byEuler(DPdata1.Euler');
    phaseDP = DPdata1.Phase;
    optDP=struct;

    % read some fields
    EBSD_fieldnames=fieldnames(EBSDdata);
    DP_EBSD_fieldnames=fieldnames(DPdata1);

    num_fields=size(EBSD_fieldnames,1);
    for n = 1: num_fields
        s=EBSDdata.(EBSD_fieldnames{n});
        if size(s,2) == 1 && size(s,1) == numel(phase)
            try
                opt.(EBSD_fieldnames{n})=s;
            catch
            end

        end
    end

    num_fields=size(DP_EBSD_fieldnames,1);
    for n = 1: num_fields
        s=DPdata1.(DP_EBSD_fieldnames{n});
        if size(s,2) == 1 && size(s,1) == numel(phase)
            try
                optDP.(DP_EBSD_fieldnames{n})=s;
            catch
            end

        end
    end

    %now put the Euler angles in the options too
    if isfield(EBSDdata,'Euler')
        opt.euler1=EBSDdata.Euler(1,:);
        opt.euler2=EBSDdata.Euler(2,:);
        opt.euler3=EBSDdata.Euler(3,:);
    end

    opt.x=optOrig.Beam_Position_X;
    opt.y=optOrig.Beam_Position_Y;

    %     opt.x=opt.X;
    %     opt.y=opt.Y;


    %now put the Euler angles in the options too
    if isfield(DPdata1,'Euler')
        optDP.euler1=DPdata1.Euler(1,:);
        optDP.euler2=DPdata1.Euler(2,:);
        optDP.euler3=DPdata1.Euler(3,:);
    end

    optDP.x=opt.x;
    optDP.y=opt.y;



    %     opt=rmfield(opt,'X');
    %     opt=rmfield(opt,'Y');

    %     opt.bc = EBSDdata.Band_Contrast;
    %     opt.bs = EBSDdata.Band_Slope;
    %     opt.bands = EBSDdata.Bands;
    %     opt.MAD = EBSDdata.Mean_Angular_Deviation;
    %     opt.quality = EBSDdata.Pattern_Quality;
    %
    % if available, add EDS data
    if exist('EDSdata','var')
        eds_names = fieldnames(EDSdata);
        for j =1 :length(eds_names)
            opt.(eds_names{j}) = EDSdata.(eds_names{j});
        end
    end

    ebsdtemp = EBSD(rot,phase,CS,opt,'unitCell',calcUnitCell([opt.x,opt.y]));
    ebsdtemp.opt.Header = EBSDheader;

    ebsdtempDP = EBSD(rotDP,phaseDP,CS,optDP,'unitCell',calcUnitCell([optDP.x,optDP.y]));
    ebsdtempDP.opt.Header = DPdata2;

    if length(EBSD_index) > 1
        ebsd{k} = ebsdtemp;
    else
        ebsd = ebsdtemp;
    end

    if length(EBSD_index) > 1
        ebsdDP{k} = ebsdtempDP;
    else
        ebsdDP = ebsdtempDP;
    end

end

% sort extra parts
ebsdDP.prop.CrossCorr=DPdata3.Cross_Correlation_Coefficient;
ebsdDP.prop.Band_Contrast=ebsdOrig.prop.Band_Contrast;
ebsdDP.prop.Band_Slope=ebsdOrig.prop.Band_Slope;
ebsdDP.prop.Pattern_Quality=ebsdOrig.prop.Pattern_Quality;
ebsdDP.prop.X=ebsdOrig.prop.X;
ebsdDP.prop.Y=ebsdOrig.prop.Y;

end

function [UCell]=Build_UCell(filename)
%Builds unit cell information from a phase file
if exist(filename) == 0
    error(['The input phase file ' filename ' cannot be accessed.']);
end

delimiter = '';
formatSpec = '%s%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);
fclose(fileID);
PhaseRaw = [dataArray{1:end-1}];

array_num=1:size(PhaseRaw,1);

ix=array_num(contains(PhaseRaw,'$name','IgnoreCase',true));
UCell.PhaseName=PhaseRaw{ix+1};

ix=array_num(contains(PhaseRaw,'$atoms','IgnoreCase',true));
atomdetails=PhaseRaw{ix+1};
semicolon_loc=strfind(atomdetails,';');
num_atoms=size(semicolon_loc,2);
semicolon_loc=[0,semicolon_loc];

for n=1:num_atoms
    atom_set=atomdetails(semicolon_loc(n)+1:semicolon_loc(n+1)-1);
    comma_loc=strfind(atom_set,',');
    
    a1_text=atom_set(1:comma_loc(1)-1);
    div_loc=strfind(a1_text,'/');
    if ~isempty(div_loc)
        a1=str2double(a1_text(1:(div_loc-1)))/str2double(a1_text((div_loc+1):end));
    else
        a1=str2double(a1_text);
    end
    
    a1_text=atom_set(comma_loc(1)+1:comma_loc(2)-1);
    div_loc=strfind(a1_text,'/');
    if isempty(div_loc)
        a2=str2double(a1_text);
    else
        a2=str2double(a1_text(1:(div_loc-1)))/str2double(a1_text((div_loc+1):end));
    end
    
    a1_text=atom_set(comma_loc(2)+1:end);
    div_loc=strfind(a1_text,'/');
    if ~isempty(div_loc)
        a3=str2double(a1_text(1:(div_loc-1)))/str2double(a1_text((div_loc+1):end));
    else
        a3=str2double(a1_text);
    end
    atom_double(n,:)=[a1,a2,a3];
    
end
UCell.atomdetails=atom_double;

ix=array_num(contains(PhaseRaw,'$symmetry','IgnoreCase',true));
UCell.crystal=str2double(PhaseRaw{ix+1});

ix=array_num(contains(PhaseRaw,'$lattice','IgnoreCase',true));
lattice_p=PhaseRaw{ix+1};
comma_loc=strfind(lattice_p,',');
UCell.a     =str2double(lattice_p(1:comma_loc(1)-1));
UCell.b     =str2double(lattice_p(comma_loc(1)+1:comma_loc(2)-1));
UCell.c     =str2double(lattice_p(comma_loc(2)+1:comma_loc(3)-1));
UCell.alpha =str2double(lattice_p(comma_loc(3)+1:comma_loc(4)-1));
UCell.beta  =str2double(lattice_p(comma_loc(4)+1:comma_loc(5)-1));
UCell.gamma =str2double(lattice_p(comma_loc(5)+1:end));

ix=array_num(contains(PhaseRaw,'$FFactor','IgnoreCase',true));
ffactors=PhaseRaw{ix+1};
comma_loc=strfind(ffactors,',');

if isempty(comma_loc)
    factor_n=str2double(ffactors);
else
    factor_n(1)=str2double(ffactors(1:comma_loc(1)-1));
    if size(comma_loc,2) > 1
        for n=2:size(comma_loc(2))-1
            factor_n(n)=str2double(ffactors(comma_loc(n-1)+1:comma_loc(n)-1));
        end
        factor_n(n+1)=str2double(ffactors(comma_loc(n)+1:end));
    end
end
UCell.ffactor=factor_n;

ix=array_num(contains(PhaseRaw,'$MaxFam','IgnoreCase',true));
UCell.max_fam=str2double(PhaseRaw{ix+1});

ix=array_num(contains(PhaseRaw,'$MaxIndex','IgnoreCase',true));
UCell.max_ind=str2double(PhaseRaw{ix+1});

ix=array_num(contains(PhaseRaw,'$EFac','IgnoreCase',true));
UCell.efac=str2double(PhaseRaw{ix+1});


%% build the final part of the cell - as per the Tutorial paper
UC.f=sqrt(1.0-( cosd(UCell.alpha)*cosd(UCell.alpha)...
    +cosd(UCell.beta)*cosd(UCell.beta)...
    +cosd(UCell.gamma)*cosd(UCell.gamma))...
    +2.0*cosd(UCell.alpha)*cosd(UCell.beta)*cosd(UCell.gamma));

%equation 2
UC.ax = UCell.a * UC.f/sind(UCell.alpha);

UC.ay = UCell.a * (cosd(UCell.gamma)-cosd(UCell.alpha)*cosd(UCell.beta))...
    /sind(UCell.alpha);

UC.az = UCell.a * cosd(UCell.beta);

%equation 3
UC.by = UCell.b * sind(UCell.alpha);
UC.bz = UCell.b * cosd(UCell.alpha);

%equation 4
UC.cz = UCell.c;

%equation 5
UCell.StructureMat=[UC.ax , 0,  0;
    UC.ay , UC.by, 0;
    UC.az , UC.bz, UC.cz];


%U.At = transpose of structure matrix as used for UVW conversions
UCell.At=transpose(UCell.StructureMat);

%U.AStar = recriprical structure matrix as used for HKL conversions
UCell.Astar=inv(UCell.StructureMat);

[UCell.num_rot,UCell.rot_m]=RotSym(UCell.crystal);

ix=array_num(contains(PhaseRaw,'$PlaneReflector','IgnoreCase',true));
if ~isempty(ix)
planedetails=PhaseRaw{ix+1};
semicolon_loc1=strfind(planedetails,';');
num_planes=size(semicolon_loc1,2);
semicolon_loc1=[0,semicolon_loc1];
plane_double=zeros(num_planes,3);
for n=1:num_planes
    plane_set=planedetails(semicolon_loc1(n)+1:semicolon_loc1(n+1)-1);
    comma_loc=strfind(plane_set,',');
    
    a1_text=plane_set(1:comma_loc(1)-1);
    div_loc=strfind(a1_text,'/');
    if ~isempty(div_loc)
        a1=str2double(a1_text(1:(div_loc-1)))/str2double(a1_text((div_loc+1):end));
    else
        a1=str2double(a1_text);
    end
    
    a1_text=plane_set(comma_loc(1)+1:comma_loc(2)-1);
    div_loc=strfind(a1_text,'/');
    if isempty(div_loc)
        a2=str2double(a1_text);
    else
        a2=str2double(a1_text(1:(div_loc-1)))/str2double(a1_text((div_loc+1):end));
    end
    
    a1_text=plane_set(comma_loc(2)+1:end);
    div_loc=strfind(a1_text,'/');
    if ~isempty(div_loc)
        a3=str2double(a1_text(1:(div_loc-1)))/str2double(a1_text((div_loc+1):end));
    else
        a3=str2double(a1_text);
    end
    plane_double(n,:)=[a1,a2,a3];
    
end
UCell.reflectors=plane_double;
UCell.max_fam=size(UCell.reflectors,1);

end

end

function [num_rot,rot_m]=RotSym(CrystalType) %generate the two symmetric matrix sets

switch CrystalType
    case 1 %cubic
        rot_m(:,:,1) = [1 0 0; 0 1 0; 0 0 1];
        rot_m(:,:,2) = [0 0 -1; 0 -1 0; -1 0 0];
        rot_m(:,:,3) = [0 0 -1; 0 1 0; 1 0 0];
        rot_m(:,:,4) = [-1 0 0; 0 1 0; 0 0 -1];
        
        rot_m(:,:,5) = [0 0 1; 0 1 0; -1 0 0];
        rot_m(:,:,6) = [1 0 0; 0 0 -1; 0 1 0];
        rot_m(:,:,7) = [1 0 0; 0 -1 0; 0 0 -1];
        rot_m(:,:,8) = [1 0 0; 0 0 1; 0 -1 0];
        
        rot_m(:,:,9) = [0 -1 0; 1 0 0; 0 0 1];
        rot_m(:,:,10) = [-1 0 0; 0 -1 0; 0 0 1];
        rot_m(:,:,11) = [0 1 0; -1 0 0; 0 0 1];
        rot_m(:,:,12) = [0 0 1; 1 0 0; 0 1 0];
        
        rot_m(:,:,13) = [0 1 0; 0 0 1; 1 0 0];
        rot_m(:,:,14) = [0 0 -1; -1 0 0; 0 1 0];
        rot_m(:,:,15) = [0 -1 0; 0 0 1; -1 0 0];
        rot_m(:,:,16) = [0 1 0; 0 0 -1; -1 0 0];
        
        rot_m(:,:,17) = [0 0 -1; 1 0 0; 0 -1 0];
        rot_m(:,:,18) = [0 0 1; -1 0 0; 0 -1 0];
        rot_m(:,:,19) = [0 -1 0; 0 0 -1; 1 0 0];
        rot_m(:,:,20) = [0 1 0; 1 0 0; 0 0 -1];
        
        rot_m(:,:,21) = [-1 0 0; 0 0 1; 0 1 0];
        rot_m(:,:,22) = [0 0 1; 0 -1 0; 1 0 0];
        rot_m(:,:,23) = [0 -1 0; -1 0 0; 0 0 -1];
        rot_m(:,:,24) = [-1 0 0; 0 0 -1; 0 -1 0];
        
    case 2
        c60=cosd(60);
        s60=sind(60);
        rot_m(:,:,1)=eye(3);
        rot_m(:,:,2)=[c60,s60,0;-s60,c60,0;0,0,1];
        rot_m(:,:,3)=[-c60,s60,0;-s60,-c60,0;0,0,1];
        rot_m(:,:,4)=[-1,0,0;0,-1,0;0,0,1];
        rot_m(:,:,5)=[-c60,-s60,0;s60,-c60,0;0,0,1];
        rot_m(:,:,6)=[c60,-s60,0;s60,c60,0;0,0,1];
        rot_m(:,:,7)=[1,0,0;0,-1,0;0,0,-1];
        rot_m(:,:,8)=[c60,s60,0;s60,-c60,0;0,0,-1];
        rot_m(:,:,9)=[-c60,s60,0;s60,c60,0;0,0,-1];
        rot_m(:,:,10)=[-1,0,0;0,1,0;0,0,-1];
        rot_m(:,:,11)=[-c60,-s60,0;-s60,c60,0;0,0,-1];
        rot_m(:,:,12)=[c60,-s60,0;-s60,-c60,0;0,0,-1];
end
num_rot=size(rot_m,3);

end
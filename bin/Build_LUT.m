function [LUTData] = Build_LUT(family_norm_list)

%link the planes together
norm_list=cat(1,family_norm_list{:});
num_planes=size(norm_list,1);

% Produce the plane normal list of triangles
plane_list=1:num_planes;
trigs=nchoosek(plane_list,3);

% Calculate the lengths for these triangles
vec1=norm_list(trigs(:,2),5:7)-norm_list(trigs(:,3),5:7);
vec2=norm_list(trigs(:,3),5:7)-norm_list(trigs(:,1),5:7);
vec3=norm_list(trigs(:,1),5:7)-norm_list(trigs(:,2),5:7);

trigs_lengths=[sqrt(dot(vec1,vec1,2)),sqrt(dot(vec2,vec2,2)),sqrt(dot(vec3,vec3,2))];
num_trigs=size(trigs,1);

% Sort the lengths and their indexes
trigs_lengths_sort=zeros(size(trigs_lengths));
trigs_sort=zeros(size(trigs));

for n=1:num_trigs
    [trigs_lengths_sort(n,:),ix]=sort(trigs_lengths(n,:));
    trigs_sort(n,:)=trigs(n,ix);
end

ratio_xy=round([trigs_lengths_sort(:,2)./trigs_lengths_sort(:,3),trigs_lengths_sort(:,1)./trigs_lengths_sort(:,3)]*1000)/1000;
perims=log10(sum(trigs_lengths_sort,2));

%remove those which are nan
perims(isnan(ratio_xy(:,1)),:)=[];
perims(isnan(ratio_xy(:,2)),:)=[];

trigs_sort(isnan(ratio_xy(:,1)),:)=[];
trigs_sort(isnan(ratio_xy(:,2)),:)=[];

ratio_xy(isnan(ratio_xy(:,1)),:) =[];
ratio_xy(isnan(ratio_xy(:,2)),:) =[];

%remove uncertainty
perims(ratio_xy(:,1) > 0.9,:) =[];
trigs_sort(ratio_xy(:,1) > 0.9,:) =[];
ratio_xy(ratio_xy(:,1) > 0.9,:) =[];

trig_data=[trigs_sort,ratio_xy,perims];

% Reduce the trig data to the unique ratio values
ratio_xy_round=round(ratio_xy*1000)/1000;

[ratio_xy_unique,ix]=unique(ratio_xy_round,'rows','first');

% Keep containers will all the potentials for each spot in the LUT
num_trig=size(ratio_xy_unique,1);

trig_group=cell(num_trig,1);

for n=1:num_trig
    trig_group{n}=trig_data(ratio_xy_round(:,1) == ratio_xy_unique(n,1) & ratio_xy_round(:,2) == ratio_xy_unique(n,2),:);
    trig_group_f{n}=unique([norm_list(trig_group{n}(:,1),1),norm_list(trig_group{n}(:,2),1),norm_list(trig_group{n}(:,3),1)],'rows');
end

LUTData.trig_group=trig_group;
LUTData.trig_group_f=trig_group_f;
LUTData.ratio_xy_unique=ratio_xy_unique;
LUTData.max_fam=size(family_norm_list,1);
LUTData.family_norm_list=family_norm_list;


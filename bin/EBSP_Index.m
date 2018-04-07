function [rotdata,band_data,err] = EBSP_Index(nhat,LUTData,thresh_trig,UCell,rot_det)
%EBSP_Index Index bands to generate a crystal orientation
%INPUTS
%nhat = input plane normals
%LUTData = Astro look up table of the potential planes
%thresh_trig = window size for searching in the LUT
%UCell = unit cell information
%rot_det = transformation to go from sample to detector
%
%OUTPUTS
%rotdata = structure with indexed information
%band_data = bands as indexed (useful for plotting)
%
%http://www.expmicromech.com

%embed in a try function so that if something goes wrong, it will continue
try
    if nargin < 5
        rot_det = eye(3);
    end
    
    num_rot2=size(nhat,1);
    % From these Y=MX+C description, calculate the normal vectors
    
    %store the useful info about the bands
    rot_list=1:num_rot2*2;
    band_trigs=nchoosek(rot_list,3);
    
    band_info=[rot_list(1:num_rot2)',nhat;[rot_list(1:num_rot2)',-nhat]];
    
    %index the bands to their families
    band_acc=zeros(num_rot2,UCell.max_fam);
    
    
    num_trig=size(band_trigs,1);
    vec1=band_info(band_trigs(:,2),2:4)-band_info(band_trigs(:,3),2:4);
    vec2=band_info(band_trigs(:,3),2:4)-band_info(band_trigs(:,1),2:4);
    vec3=band_info(band_trigs(:,1),2:4)-band_info(band_trigs(:,2),2:4);
    
    band_trigs_lengths=[sqrt(dot(vec1,vec1,2)),sqrt(dot(vec2,vec2,2)),sqrt(dot(vec3,vec3,2))];
    
    %sort the lengths
    band_trigs_lengths_sort=zeros(size(band_trigs_lengths));
    band_trigs_sort=zeros(size(band_trigs));
    
    for n=1:num_trig
        [band_trigs_lengths_sort(n,:),ix]=sort(band_trigs_lengths(n,:));
        band_trigs_sort(n,:)=band_trigs(n,ix);
    end
    
    band_ratio_xy=round([band_trigs_lengths_sort(:,2)./band_trigs_lengths_sort(:,3),band_trigs_lengths_sort(:,1)./band_trigs_lengths_sort(:,3)]*1000)/1000;
    band_perims=log10(sum(band_trigs_lengths_sort,2));
    
    band_trig_data=[band_trigs_sort,band_ratio_xy,band_perims];
    
    [band_ratio_xy_unique,ix]=unique(band_ratio_xy,'rows','first');
    
    band_ratio_xy_unique(band_ratio_xy_unique(:,1) > 0.9,:) =[];
    band_ratio_xy_unique(isnan(band_ratio_xy_unique(:,1)),:) =[];
    band_ratio_xy_unique(isnan(band_ratio_xy_unique(:,2)),:) =[];
    band_ratio_xy_unique(round (band_ratio_xy_unique(:,1) * UCell.efac) == round (band_ratio_xy_unique(:,2) * UCell.efac),:) = [];
    
    num_trig=size(band_ratio_xy_unique,1);
    
    %from the unique ratios, work out the groupings that correspond to them
    
    band_trig_group=cell(num_trig,1);
    num_band_trig=zeros(num_trig,1);
    
    num_xyband=size(band_ratio_xy_unique,1);
    
    band_trig_data(band_trig_data(:,1)>num_rot2,1)=band_trig_data(band_trig_data(:,1)>num_rot2,1)-num_rot2;
    band_trig_data(band_trig_data(:,2)>num_rot2,2)=band_trig_data(band_trig_data(:,2)>num_rot2,2)-num_rot2;
    band_trig_data(band_trig_data(:,3)>num_rot2,3)=band_trig_data(band_trig_data(:,3)>num_rot2,3)-num_rot2;
    
    band_trig_data=sortrows(unique(band_trig_data,'rows','first'),[4,5,6]);
    
    for n=1:num_xyband
        
        band_trig_group{n}=band_trig_data(band_trig_data(:,4) == band_ratio_xy_unique(n,1) & band_trig_data(:,5) == band_ratio_xy_unique(n,2),:);
        num_band_trig(n)=size(band_trig_group{n},1);
        
        %now calculate the candidates
        %reduce the search space
        reduced_index=(LUTData.ratio_xy_unique(:,2)> band_ratio_xy_unique(n,2) - thresh_trig) & (LUTData.ratio_xy_unique(:,2)< band_ratio_xy_unique(n,2) + thresh_trig);
        reduced_xylut=LUTData.ratio_xy_unique(reduced_index,:);
        reduced_trig_group_f=LUTData.trig_group_f(reduced_index);
        
        %calculate the distances
        reduced_distance=sqrt((reduced_xylut(:,1)-band_ratio_xy_unique(n,1)).^2+(reduced_xylut(:,1)-band_ratio_xy_unique(n,1)).^2);
        reduced2_index=(reduced_distance < thresh_trig);
        candidate_trig_cell=reduced_trig_group_f(reduced2_index);
        candidate_trig=cat(1,candidate_trig_cell{:});
        
        %extract the family accumulator points for each band 'parent' (i.e. band = p)
        if isempty(candidate_trig) == 0
            for s=1:num_band_trig(n)
                for t=1:3
                    ax=band_trig_group{n}(s,t);
                    band_acc(ax,candidate_trig(:,t))=band_acc(ax,candidate_trig(:,t))+1;
                end
            end
        end
        
    end
    
    
    %calculate the trust factor and most likely family
    [band_max,band_id]=max(band_acc,[],2);
    
    
    band_id_poss=nchoosek(1:size(band_id,1),2);
    num_choice=size(band_id_poss,1);
    
    band_oknum=zeros(size(band_id_poss,1),1);
    band_fid_ang=zeros(num_choice,size(band_id,1));
    band_fid=zeros(num_choice,size(band_id,1));
    band_ref=zeros(3,size(band_id,1),num_choice);
    band_ok_group=zeros(size(band_id,1),num_choice);
    
    for p=1:num_choice
        band_id_start=band_id_poss(p,:);
        
        checkdot_temp=acosd(nhat(band_id_start(1),:)*nhat(band_id_start(2),:).');
        if checkdot_temp > 5 %not colinear
            %calculate their interplanar angle
            ang_12=nhat(band_id_start(1),:)*nhat(band_id_start(2),:).';
            
            %select the first band from the first family
            band1_norm_ref=LUTData.family_norm_list{band_id(band_id_start(1))}(1,5:7);
            %create a list of bands from the second family
            band2_norm_ref_poss=LUTData.family_norm_list{band_id(band_id_start(2))}(:,5:7);
            %calculate possible interplanar angles
            
            %check the angluar different between the measured pair and the list
            [~,band2_id]=min(abs(band1_norm_ref*band2_norm_ref_poss'-ang_12));
            
            %pick the closest match
            band2_norm_ref=band2_norm_ref_poss(band2_id,:);
            
            %We now have two bands which are from a suitable set of families of planes
            %Now solve for the rough orientation matrix
            
            %Calculate the normal to band 1 and band 2
            %from the reference orientation
            %         b12_nref=cross(band1_norm_ref,band2_norm_ref);
            
            %cross product - quicker to calculate manually
            
            b12_nref=[band1_norm_ref(2)*band2_norm_ref(3)-band1_norm_ref(3)*band2_norm_ref(2),...
                band1_norm_ref(3)*band2_norm_ref(1)-band1_norm_ref(1)*band2_norm_ref(3),...
                band1_norm_ref(1)*band2_norm_ref(2)-band1_norm_ref(2)*band2_norm_ref(1)];
            
            %normalised length
            b12_nref=b12_nref(:)/(norm(b12_nref));
            
            %from the test orientation
            
            %cross product - quicker to calculate manually
            b12_ntest=[nhat(band_id_start(1),2)*nhat(band_id_start(2),3)-nhat(band_id_start(1),3)*nhat(band_id_start(2),2),...
                nhat(band_id_start(1),3)*nhat(band_id_start(2),1)-nhat(band_id_start(1),1)*nhat(band_id_start(2),3),...
                nhat(band_id_start(1),1)*nhat(band_id_start(2),2)-nhat(band_id_start(1),2)*nhat(band_id_start(2),1)];
            
            b12_ntest=b12_ntest(:)/norm(b12_ntest);
            
            %construct the three vector matricies
            p_ref=[band1_norm_ref(:),band2_norm_ref(:),b12_nref];
            p_test=[nhat(band_id_start(1),:).',nhat(band_id_start(2),:).',b12_ntest];
            
            %calculate the rotation matrix (rough - only from these two bands
            %and not an orthogonal transformation!
            mat_rot=p_ref/p_test;
            
            %now that the transformation has been determined, we can 'index' all the other bands...
            
            bands_inref=mat_rot*nhat.';
            
            for n=1:size(bands_inref,2)
                ref_bandlist=LUTData.family_norm_list{band_id(n)}(:,5:7).';
                [band_fid_ang(p,n),band_fid(p,n)]=max(ref_bandlist'*bands_inref(:,n));
                band_ref(:,n,p)=ref_bandlist(:,band_fid(p,n));
            end
            
            band_keep=find(band_fid_ang(p,:) > 0.996194698091746); %cosd3
            band_ok=false(1,size(bands_inref,2));
            band_ok(band_keep)=true; %#ok<FNDSB>
            band_oknum(p)=sum(band_ok);
            band_ok_group(:,p)=band_ok;
        else
            band_oknum(p)=0;
        end
    end
    
    %pick the solution that has the most bands matching
    [band_maxnum,band_maxok_ind]=max(band_oknum);
    
    %from these indexed bands, calculate the orientation matrix
    %(using all the bands)
    
    band_data=[band_id,band_max,nhat,band_ref(:,:,band_maxok_ind).',band_fid_ang(band_maxok_ind,:).',band_fid(band_maxok_ind,:).',band_ok_group(:,band_maxok_ind)];
    %[1 = band family, 2 - max points, 3-5 = nhat rotated, 6-8
    ok_bands=logical(band_ok_group(:,band_maxok_ind));
    
    b_ref=squeeze(band_ref(:,ok_bands,band_maxok_ind)).';
    b_test=nhat(ok_bands,:);
    %calc a transformation matrix close to the rotation matrix
    if size(b_test,1) > 3
        b_test2=b_test-repmat(mean(b_test,1),[size(b_test,1) 1]);
        b_ref2=b_ref-repmat(mean(b_ref,1),[size(b_ref,1) 1]);
        %use the Kabash algorithm
        co_var=b_test2.'*b_ref2;
        [a,~,c]=svd(co_var);
        d=sign(det(a'*c));
        rotdata.detector=c*[1,0,0;0,1,0;0,0,d]*a.';
        rotdata.sample=rotdata.detector*rot_det';
        rotdata.eang=conv_G_to_EA(rotdata.sample);
        %back calc
        rotdata.bfinal=b_ref*rotdata.detector;
        %calc the error
        rotdata.error=mean(acos(dot(rotdata.bfinal,b_test,2)));
        rotdata.maxok=band_maxnum;
        
    else
        error('Not enough bands matching')
    end
catch err
    rotdata.eang=[0 0 0];
    rotdata.detector=eye(3);
    rotdata.sample=eye(3);
    %calculat the different between the regularised version and the transform
    rotdata.error=nan;
    rotdata.maxok=nan;
    rotdata.bfinal=nan;
    band_data=[0];
end



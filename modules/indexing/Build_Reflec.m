function [ Family_List ] = Build_Reflec( UCell )
%GEN_REFLEC3 Summary of this function goes here
%   Detailed explanation goes here

randv=[0.434020560016200   0.577508235281481   0.691455270907160];

if isfield(UCell,'reflectors') == 0
    index_list=0:UCell.max_ind;
    [h_grid,k_grid,l_grid]=meshgrid(index_list,index_list,index_list);
    HKL_list=[h_grid(:),k_grid(:),l_grid(:)];
    HKL_list(1,:)=[]; %remove [0 0 0];
    
    HKL_list_cart=HKL_list*UCell.Astar; %check that this is the right equation...
    HKL_list_len=round(1./sqrt(dot(HKL_list_cart,HKL_list_cart,2))*100000)/100000;
    
    %reduce remove higher order
    [HKL_len_unique,ix]=unique(HKL_list_len);
    ix=flipud(ix);
    HKL_list_ix=HKL_list(ix,:);
%     HKL_list_cart_ix=HKL_list_cart(ix,:);
    
    %calculate the structure factors
    sfactor=zeros(size(ix,1),1);
    for n=1:size(ix,1)
        sfactora=norm(sum(UCell.ffactor(:).*exp(2*pi*1i*(UCell.atomdetails(:,1)*HKL_list_ix(n,1)+UCell.atomdetails(:,2)*HKL_list_ix(n,2)+UCell.atomdetails(:,3)*HKL_list_ix(n,3)))));
        sfactor(n)=round((sfactora)*100000)/100000;
    end
    
    %remove the invisible planes
    ix1=ix(sfactor>0);
    HKL_list_ix=HKL_list(ix1,:);
    HKL_list_cart_ix=HKL_list_cart(ix1,:);
    HKL_list_len_ix=HKL_list_len(ix1);
    
    %dump the higher order reflectors
    angle_re=dot(repmat(randv,size(HKL_list_cart_ix,1),1),repmat(HKL_list_len_ix,1,3).*HKL_list_cart_ix,2);
    [~,ix2]=unique(round(angle_re*100000)/100000);
    ix3=ix1(ix2);
    
%     HKL_list_ix=HKL_list(ix3,:);
%     HKL_list_cart_ix=HKL_list_cart(ix3,:);
    HKL_list_len_ix=HKL_list_len(ix3);
    
    %sort for plane spacing
    [~,ix4]=sort(HKL_list_len_ix,'descend');
    ix4=ix3(ix4);
    
%     HKL_list_ix=HKL_list(ix4,:);
%     HKL_list_cart_ix=HKL_list_cart(ix4,:);
%     HKL_list_len_ix=HKL_list_len(ix4);
    
    %need to apply symmetry within each family % find unique terms
    ix5=ix4(1:UCell.max_fam);
%     HKL_list_ix=HKL_list(ix5,:);
    HKL_list_cart_ix=HKL_list_cart(ix5,:);
%     HKL_list_len_ix=HKL_list_len(ix5);
    
    Family_List=cell(UCell.max_fam,1);
else
    UCell.max_fam=size(UCell.reflectors,1);
    Family_List=cell(UCell.max_fam,1);
    HKL_list_cart_ix=UCell.reflectors*UCell.Astar; %check that this is the right equation
end

for n=1:UCell.max_fam
    HKL_list_cart_r=zeros(UCell.num_rot,3);
    for s=1:UCell.num_rot
        HKL_list_cart_r(s,:)=HKL_list_cart_ix(n,:)*UCell.rot_m(:,:,s);
    end
    [~,ik]=unique(round(HKL_list_cart_r*10000)/10000,'rows');
    HKL_list_cart_n=HKL_list_cart_r(ik,:)./repmat(sqrt(dot(HKL_list_cart_r(ik,:),HKL_list_cart_r(ik,:),2)),1,3);
    HKL_list_f=HKL_list_cart_r(ik,:)*inv(UCell.Astar);
    num_f=size(ik,1);
    %       Bandset=[repmat([n,HKL_list_len_ix(n),HKL_list_ix(n,:)],[num_f,1]),HKL_list_cart_f,HKL_list_cart_n,HKL_list_f];
    Bandset=[repmat(n,size(HKL_list_f,1),1), HKL_list_f,      HKL_list_cart_n];
    Family_List{n}=Bandset;
end
%Bandset= [family number, HKL_len, HKL_cart,rotated cart,rotated HKL (rounded)]

end


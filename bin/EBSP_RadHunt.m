function [ Peak_Centre,Peak_Quality,Peak_Set,R_EBSP,R_Edge,R_RHO,R_theta ] = EBSP_RadHunt( EBSP,Settings_Rad )
%EBSP_RadHunt Radon transform the EBSP & find peaks
%[ Peak_Centre,Peak_Quality,Peak_Set,R_EBSP,R_Edge,R_RHO,R_Theta ] = EBSP_RadHunt( EBSP,Settings_Rad )
%
%INPUTS
%EBSP = input EBSP pattern (greyscale, double)
%Settings_Rad = Radon Transform Settings (as structure)
%
%OUTPUTS
%Peak_Centre = [theta,rho,n] or peak centres
%Peak_Quality =[IQ,Slope]
%Peak_Set = [theta_edge1,rho_edge1,height1,theta_edge2,rho_edge2,heigh2]
%R_EBSP = transformed EBSP
%R_Edge = Sobel filtered EBSP
%R_RHO = rhos
%R_theta = thetas

%% Versioning
%v1 - TBB 14/04/2017

%%

%transform
R_theta=Settings_Rad.theta_range(1):Settings_Rad.theta_range(3):Settings_Rad.theta_range(2);
[R_EBSP,R_RHO]=radon(EBSP,R_theta);
[R_EBSP_b,R_RHO]=radon(EBSP*0+1,R_theta); %creates the zero mask
R_EBSP=R_EBSP./R_EBSP_b;
%remove the areas outside the mask
R_EBSP(R_EBSP==0)=NaN;
R_EBSP(R_EBSP_b==0)=NaN;
R_RHO_num=numel(R_RHO);

Settings_Rad.rho_search_pix=Settings_Rad.rho_search_per*R_RHO_num;

%set up the search space
theta_box=-Settings_Rad.theta_search_pix:Settings_Rad.theta_search_pix;
theta_box2=-3*Settings_Rad.theta_search_pix:3*Settings_Rad.theta_search_pix;
rho_box=round(-Settings_Rad.rho_search_pix:0);
Settings_Rad.rho_search_pix_plus=floor(Settings_Rad.rho_search_pix/2);
rho_box2=round(-Settings_Rad.rho_search_pix:Settings_Rad.rho_search_pix_plus);


%conv the radon
edgefilter=[-1,-2,-1; 0 0 0; 1 2 1];
R_Edge=conv2(R_EBSP,edgefilter,'same');

Peak_X1=zeros(Settings_Rad.num_peak,1);
Peak_Y1=zeros(Settings_Rad.num_peak,1);
Peak_X2=zeros(Settings_Rad.num_peak,1);
Peak_Y2=zeros(Settings_Rad.num_peak,1);
Peak_H1=zeros(Settings_Rad.num_peak,1);
Peak_H2=Peak_H1;

%downscale factor for peaks
% dfac=1E-9;
dfac=0;

Pat_RadSob=R_Edge;
x_size=size(R_Edge,2);
y_size=size(R_Edge,1);

for n=1:Settings_Rad.num_peak
    %find the max value
    [Yvals,Ypos3]=nanmax(Pat_RadSob);
    [Peak_H1(n),Peak_X1(n)]=nanmax(Yvals);
    Peak_Y1(n)=Ypos3(Peak_X1(n));
    
    x_box =Peak_X1(n)+theta_box;
    x_box2=Peak_X1(n)+theta_box2;
    y_box =Peak_Y1(n)+rho_box;
    y_box2=Peak_Y1(n)+rho_box2;
    x_box(x_box<=0 | x_box>x_size) = [];
    x_box2(x_box2<=0 | x_box2>x_size) = [];
    y_box(y_box<=0 | y_box>y_size) = [];
    y_box2(y_box2<=0 | y_box2>y_size) = [];
    
    %carve out the box
    SubWindow=Pat_RadSob(y_box,x_box);
    
    if size(x_box,2) == size(theta_box,2) && size(y_box,2) == size(rho_box,2) && size(y_box2,2) == size(rho_box2,2)
        [Yvals,Ypos3]=nanmin(SubWindow);
        [Peak_H2(n),Peak_X]=nanmin(Yvals);
        Peak_Y=Ypos3(Peak_X);
        
        Peak_Y2(n)=Peak_Y1(n)+rho_box(Peak_Y);
        Peak_X2(n)=Peak_X1(n)+theta_box(Peak_X);
    else %box is for a near end of screen peak
        Peak_H1(n)=0;
        Peak_Y1(n)=1;
        Peak_X1(n)=1;
        Peak_H2(n)=0;
        Peak_Y2(n)=1;
        Peak_X2(n)=1;
    end
        
    %downweight the area for the next peak search
    Pat_RadSob(y_box2,x_box)=Pat_RadSob(y_box2,x_box)*dfac;
end

%eR_XPort the peaks found
Peak_Set=[R_theta(Peak_X1).',R_RHO(Peak_Y1),Peak_H1,R_theta(Peak_X2).',R_RHO(Peak_Y2),Peak_H2];

%remove the duplicates
[~,Peak_Unique,~]=unique(round(Peak_Set(:,1:2)*1000)/1000,'rows');
Peak_Set=Peak_Set(Peak_Unique,:);

Peak_Centre=[Peak_Set(:,1)+Peak_Set(:,4),Peak_Set(:,2)+Peak_Set(:,5)]./2;
Peak_Width=Peak_Set(:,2)-Peak_Set(:,5);

Peak_Ok=find(Peak_Centre(:,1) >= -2 & Peak_Centre(:,1) <= 180 & Peak_Width > Settings_Rad.min_peak_width*R_RHO_num);

Peak_Centre=Peak_Centre(Peak_Ok,:); %remove the repeats
Peak_Set=Peak_Set(Peak_Ok,:); %remove the repeats

%need the index for each peak centre to get the peak height
%Peak_Centre(:,1) %thetas
%Peak_Centre(:,2) %rhos

%extract the peak centres
[~,itheta]=min(abs(repmat(Peak_Centre(:,1),1,numel(R_theta))-repmat(R_theta,numel(Peak_Centre(:,1)),1)),[],2);
[~,irho]  =min(abs(repmat(Peak_Centre(:,2),1,numel(R_RHO))-repmat(R_RHO',numel(Peak_Centre(:,2)),1)),[],2);
ilist=sub2ind(size(R_EBSP),irho,itheta);
pheight=R_EBSP(ilist);

%pick peaks of highest BC
edge_height= Peak_Set(:,3)-Peak_Set(:,6);
[~,ix]=sort(edge_height,'descend');

max_peaks=min([Settings_Rad.max_peaks numel(ix)]);

ix=ix(1:max_peaks);

%select only sharp bands
edge_height=edge_height(ix);
pheight=pheight(ix);
Peak_Centre=Peak_Centre(ix,:);
Peak_Set=Peak_Set(ix,:);


%Image quality metric (band sharpness, adapted from Kunze 1993)
edges=(nansum(edge_height)/numel(edge_height))/nanstd(R_Edge(:));
Peak_Quality=[mean(pheight) edges];


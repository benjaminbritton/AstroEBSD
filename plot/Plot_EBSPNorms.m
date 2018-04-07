function Plot_EBSPNorms(EBSD_Pattern_norm,EBSD_Geometry_in,nhat_gnom,s1 )
%PLOT_EBSPNORMS Summary of this function goes here
%   Detailed explanation goes here

x2_gnom=zeros(size(nhat_gnom,1),2);
y2_gnom=zeros(size(nhat_gnom,1),2);

for n=1:size(nhat_gnom,1)
x2_gnom(n,:)=[EBSD_Geometry_in.x_screen(1) EBSD_Geometry_in.x_screen(end)];
y2_gnom(n,:)=(-x2_gnom(n,:).*nhat_gnom(n,1)-nhat_gnom(n,3))./nhat_gnom(n,2);
end

imagesc(EBSD_Geometry_in.x_screen,EBSD_Geometry_in.y_screen,EBSD_Pattern_norm,'Parent',s1);
    
for n=1:size(nhat_gnom)
    hold on;
%     plot(x1_gnom(n,:),y1_gnom(n,:),'r');
    plot(x2_gnom(n,:),y2_gnom(n,:),':','color','r','LineWidth',2,'Parent',s1);
% plot(x2_gnom(n,:),y2_gnom(n,:),'.','color','w','LineWidth',2,'Parent',s1);
end
% 
% xlabel(s1,'X / Z');
% ylabel(s1,'Y / Z');

% colormap('gray')
% axis equal; axis xy;
% ylim([EBSD_Geometry_in.y_screen(1) EBSD_Geometry_in.y_screen(end)]);
% xlim([EBSD_Geometry_in.x_screen(1) EBSD_Geometry_in.x_screen(end)]);

end


function [ f1 ] = Plot_EBSPAnnotated( EBSD_Pattern_in,EBSD_Geometry,nhat_gnom,RotMat,Crystal_UCell,Family_List,s1)
%PLOT_EBSPANNOTATED Summary of this function goes here
%   Detailed explanation goes here
cset=[1,0,0;
    0,1,0;
    0,0,1;
    0,1,1;
    1,0,1;
    1,1,0];
cset=[cset;cset*0.5];

U.K         =   Crystal_UCell.At*RotMat;
U.Kstar     =   Crystal_UCell.Astar*RotMat;


f1(1)=imagesc(EBSD_Geometry.x_screen,EBSD_Geometry.y_screen,EBSD_Pattern_in,'Parent',s1);

axis(s1,'equal','xy');
hold(s1,'on')
%label the axes
xlabel(s1,'X / Z');
ylabel(s1,'Y / Z');

colormap('gray')

%plot the PC
f1(2)=scatter(0,0,50,'wo','filled','Parent',s1);
f1(3)=scatter(0,0,50,'rx','Parent',s1);




% plot the bands from the rotated solution

if ~isempty(Family_List)
    for n=1:size(Family_List,1)
        HKLs=Family_List{n}(:,5:7)*RotMat;
        
        for p=1:size(HKLs,1)
            nhat=HKLs(p,:);
            nhat=nhat/norm(nhat);
            
            if abs(nhat(2)) > 0.5
                xlin=linspace(EBSD_Geometry.x_screen(1),EBSD_Geometry.x_screen(end),100);
                ylin=(-xlin.*nhat(1)-nhat(3))./nhat(2);
                
                [ix]=find(ylin > EBSD_Geometry.y_screen(1) & ylin < EBSD_Geometry.y_screen(end));
                ylin2=ylin(ix);
                xlin2=xlin(ix);
                
                if ~isempty(ix)
                    f1(end+1)=plot(xlin2,ylin2,'color',cset(n,:),'LineWidth',2,'parent',s1);
                end
            else
                
                ylin=linspace(EBSD_Geometry.y_screen(1),EBSD_Geometry.y_screen(end),100);
                xlin=(-ylin.*nhat(2)-nhat(3))./nhat(1);
                
                [ix]=find(xlin > EBSD_Geometry.x_screen(1) & xlin < EBSD_Geometry.x_screen(end));
                ylin2=ylin(ix);
                xlin2=xlin(ix);
                
                if ~isempty(ix)
                    f1(end+1)=plot(xlin2,ylin2,'color',cset(n,:),'LineWidth',2,'Parent',s1);
                end
                
            end
            
        end
    end
end




%plot the detected bands
if ~isempty(nhat_gnom)
    for n=1:size(nhat_gnom,1)
        
        nhat=nhat_gnom(n,:);
        nhat=nhat/norm(nhat);
        
        if abs(nhat(2)) > 0.5
            xlin=linspace(EBSD_Geometry.x_screen(1),EBSD_Geometry.x_screen(end),100);
            ylin=(-xlin.*nhat(1)-nhat(3))./nhat(2);
            
            [ix]=find(ylin > EBSD_Geometry.y_screen(1) & ylin < EBSD_Geometry.y_screen(end));
            ylin2=ylin(ix);
            xlin2=xlin(ix);
            
            if ~isempty(ix)
                f1(end+1)=plot(xlin2,ylin2,'.','color','w','LineWidth',2,'Parent',s1);
            end
        else
            
            ylin=linspace(EBSD_Geometry.y_screen(1),EBSD_Geometry.y_screen(end),100);
            xlin=(-ylin.*nhat(2)-nhat(3))./nhat(1);
            
            [ix]=find(xlin > EBSD_Geometry.x_screen(1) & xlin < EBSD_Geometry.x_screen(end));
            ylin2=ylin(ix);
            xlin2=xlin(ix);
            
            if ~isempty(ix)
                f1(end+1)=plot(xlin2,ylin2,'.','color','w','LineWidth',2,'Parent',s1);
            end
            
        end
        
    end
end

ylim(s1,[EBSD_Geometry.y_screen(1) EBSD_Geometry.y_screen(end)]);
xlim(s1,[EBSD_Geometry.x_screen(1) EBSD_Geometry.x_screen(end)]);
end


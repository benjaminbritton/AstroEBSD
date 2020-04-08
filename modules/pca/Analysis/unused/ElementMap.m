%This script defines some peak heights found in Co alloys and grabs their
%intensity for each RC-spectrum.
%Then plots some chemical maps.

%TPM 07/06/19

% Elemental peak maps
w(:,1)=[40,63];
w(:,2)=[66,91];%C
w(:,3)=[100,113];%Cr
w(:,4)=[125,136]; %two peaks: Co, Ni
w(:,5)=[191,204]; %Al
w(:,6)=[220,233]; % Ta/W
w(:,7)=[240,262]; %Zr
w(:,8)=[272,287]; %Mo
w(:,9)=[579,598]; %Cr
w(:,10)=[634,670]; %Cr
w(:,11)=[727,754]; %Co
w(:,12)=[784,810]; %Ni
w(:,13)=[862,890]; %Ta/W
w(:,14)=[1003,1024]; %Ta/W

EDSWindowing.elements={"-","C","Cr","Co/Ni","Al","Ta/W","Zr","Mo","Cr","Cr","Co","Ni","Ta/W","Ta/W"};


EDSWindowing.w=w;
clear w;

for b=1:total_ret_comps
    spectrum=tile.spec_reshaped(:,b);
    spectrum=spectrum./std(spectrum);
    
    for a=1:14
    ind1=EDSWindowing.w(1,a);
    ind2=EDSWindowing.w(2,a);
    [pk,loc]=findpeaks(spectrum(ind1:ind2),'SortStr','descend','NPeaks',1);
    
    if isempty(pk)==1
        continue
    end
    
    EDSWindowing.p(a,b)=pk;
    EDSWindowing.l(a,b)=loc;
    
    EDSWindowing.l(a,b)=EDSWindowing.l(a,b)+ind1;
    end
    
    clear spectrum ind1 ind2 pk loc
end
clear a b

for a=1:14
EDSMap=zeros(size(tile.map_reshaped));
for b=1:total_ret_comps
    locs=find(tile.map_reshaped==b);
    EDSMap(locs)=EDSWindowing.p(a,b);
    clear locs
end

EDSWindowing.Maps(:,:,a)=EDSMap;
clear EDSMap
end
clear a b

cmap=cbrewer('div','PuOr',100);

figure('Position', get(groot,'ScreenSize'))
[hsub,wsub]=tight_subplot(2,4,0.05,0.05);

axes(hsub(1))
imagesc(EDSWindowing.Maps(:,:,2))
pbaspect([1,1,1])
colormap(gca,cmap)
title('Carbon')
axes(hsub(2))
imagesc(EDSWindowing.Maps(:,:,5))
pbaspect([1,1,1])
colormap(gca,cmap)
title('Aluminium')
axes(hsub(3))
imagesc(EDSWindowing.Maps(:,:,7))
pbaspect([1,1,1])
colormap(gca,cmap)
title('Zirconium')
axes(hsub(4))
imagesc(EDSWindowing.qualmap)
pbaspect([1,1,1])
colormap(gca,cmap)
title('Similarity index')
axes(hsub(5))
imagesc(EDSWindowing.Maps(:,:,8))
pbaspect([1,1,1])
colormap(gca,cmap)
title('Molybdenum')
axes(hsub(6))
imagesc(EDSWindowing.Maps(:,:,11))
pbaspect([1,1,1])
colormap(gca,cmap)
title('Cobalt')
axes(hsub(7))
imagesc(EDSWindowing.Maps(:,:,12))
pbaspect([1,1,1])
colormap(gca,cmap)
title('Nickel')
axes(hsub(8))
imagesc(EDSWindowing.Maps(:,:,10))
pbaspect([1,1,1])
colormap(gca,cmap)
title('Chromium')

set(hsub(1:8),'XTickLabel',''); set(hsub,'YTickLabel','');
%colormap(gca,cmap)

if printing==1
    cd(InputUser.ResultsDir)
    print(gcf,['ElementMap'],'-dpng','-r300');
end



clear hsub wsub cmap

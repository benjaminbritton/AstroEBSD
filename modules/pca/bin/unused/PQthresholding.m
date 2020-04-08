pq=Plots.RadonCropped(:,:,1);
[thresh,metric]=multithresh(pq,1);
pqseg=imquantize(pq, thresh);
imagesc(pqseg)

props=regionprops(pqseg);
%%
imagesc(pq)

%%
pqbw=imbinarize(pq);
imagesc(pqbw)

%%
gmag=imgradient(pq,'sobel');
imagesc(gmag)
gmagbc=imcomplement(imbinarize(gmag));
imagesc(gmagbc)
%%
bw2 = ~bwareaopen(~gmagbc, 10);
D = -bwdist(~gmagbc);
Ld = watershed(D);
imshow(label2rgb(Ld))

%%
I2 = imcomplement(gmag.*255);
I3 = imhmin(I2,40); %20 is the height threshold for suppressing shallow minima
L = watershed(I3);
imagesc(L)

%%
I=pq;

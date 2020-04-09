function [top]=Tile_Extend(topchunk,Data_InputMap)

top=zeros(length(topchunk{1}),length(topchunk{2}));

%load the extra patterns required
for i=1:length(topchunk{1}) %loop over topchunk rows
    for j=1:length(topchunk{2})
        row=topchunk{1}(i);
        col=topchunk{2}(j);
        
        if row>0 & col>0 & row<Data_InputMap.ypts+1 & col<Data_InputMap.xpts+1
            top(i,j)=Data_InputMap.PMap(topchunk{1}(i),topchunk{2}(j));

        elseif row>0 & col>0 & row<Data_InputMap.ypts+1 & col>=Data_InputMap.xpts+1
            top(i,j)=Data_InputMap.PMap(topchunk{1}(i),Data_InputMap.xpts);

        elseif row>0 & col>0 & row>=Data_InputMap.ypts+1 & col<Data_InputMap.xpts+1
            top(i,j)=Data_InputMap.PMap(Data_InputMap.ypts,topchunk{2}(j));

        elseif row>0 & col<=0 & row<Data_InputMap.ypts+1 & col<Data_InputMap.xpts+1
            top(i,j)=Data_InputMap.PMap(topchunk{1}(i),1);

        elseif row<=0 & col>0 & row<Data_InputMap.ypts+1 & col<Data_InputMap.xpts+1
            top(i,j)=Data_InputMap.PMap(1,topchunk{2}(j));
            
        elseif row<=0 & col<=0
            top(i,j)=Data_InputMap.PMap(1,1);
            
        elseif row>=Data_InputMap.ypts+1 & col>=Data_InputMap.xpts+1
            top(i,j)=Data_InputMap.PMap(Data_InputMap.ypts,Data_InputMap.xpts);
            
        elseif row<=0 & col>=Data_InputMap.xpts+1
            top(i,j)=Data_InputMap.PMap(1,Data_InputMap.xpts);
            
        elseif row>=Data_InputMap.ypts+1 & col<=0
            top(i,j)=Data_InputMap.PMap(Data_InputMap.ypts,1);
            
            
        end
    end
end
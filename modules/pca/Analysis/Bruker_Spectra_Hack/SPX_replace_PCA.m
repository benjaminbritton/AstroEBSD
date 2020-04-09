% - T P McAuliffe 20/02/2020 for AstroEBSD v2

% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
% OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
% NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
% BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN 
% ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
% CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

% folder is the location of the original spx file
% file_original is the name of the original spx file
% spec is self explanatory
% folder 2 is destination folder for spx files

function []=SPX_replace_PCA(folder,file_original,spec,folder2)

%% Replace spx spectrum with a PCA one...
% TPM 9/1/18: this code replaces the <Channels> line of a .spx file with a
% user-inputted 'spectrum' variable. Modifies an existing .spx file.

cd(folder)
file=[file_original(1:end-4),'_modified','.spx'];
copyfile(file_original, file)

chanlength=2039;

cd(folder2)
%test spectrum is random integers...
%spectrum=randi(500,1,length);

%%
[a,b]=size(spec);

%%
for c=1:b

spectrum=spec(1:chanlength,c);
sumspectrum=spec(1:chanlength,c);
meanspectrum=spec(1:chanlength,c);
%pcaspectrum=Summing.pcaspec(1:chanlength,c);
pcaspectrum=spec(1:chanlength,c);

%%%%%%% Spectrum %%%%%%%
reinsert=sprintf('%d,',spectrum);
reinsert2=['<Channels>',reinsert,'</Channels>'];
% Read txt into cell A
fid = fopen(fullfile(folder,file),'r');
i = 1;
tline = fgetl(fid);
A{i} = tline;
while ischar(tline)
    i = i+1;
    tline = fgetl(fid);
    A{i} = tline;
end
fclose(fid);
% Change cell at position 105 (<Channels>)
A{105} = sprintf('%s',reinsert2);
% Write cell A into txt
fid = fopen(file, 'w');
for i = 1:numel(A)
    if A{i+1} == -1
        fprintf(fid,'%s', A{i});
        break
    else
        fprintf(fid,'%s\n', A{i});
    end
end
fclose('all');
movefile(file,['spectrum_',num2str(c),'.spx']);
clear ans fid i reinsert reinsert2 tline

end
% %%%%%%% Sumspectrum %%%%%%%
% reinsert=sprintf('%d,',sumspectrum);
% reinsert2=['<Channels>',reinsert,'</Channels>'];
% % Read txt into cell A
% fid = fopen([folder,file],'r');
% i = 1;
% tline = fgetl(fid);
% A{i} = tline;
% while ischar(tline)
%     i = i+1;
%     tline = fgetl(fid);
%     A{i} = tline;
% end
% fclose(fid);
% % Change cell at position 105 (<Channels>)
% A{105} = sprintf('%s',reinsert2);
% % Write cell A into txt
% fid = fopen(file, 'w');
% for i = 1:numel(A)
%     if A{i+1} == -1
%         fprintf(fid,'%s', A{i});
%         break
%     else
%         fprintf(fid,'%s\n', A{i});
%     end
% end
% mkdir('sum')
% cd('sum')
% movefile([folder3,'/',file],['sumspectrum_',num2str(c),'.spx']);
% clear ans fid i reinsert reinsert2 tline
% cd(folder3)
% 
% %%%%%%% Meanspectrum %%%%%%%
% reinsert=sprintf('%d,',meanspectrum);
% reinsert2=['<Channels>',reinsert,'</Channels>'];
% % Read txt into cell A
% fid = fopen([folder,file],'r');
% i = 1;
% tline = fgetl(fid);
% A{i} = tline;
% while ischar(tline)
%     i = i+1;
%     tline = fgetl(fid);
%     A{i} = tline;
% end
% fclose(fid);
% % Change cell at position 105 (<Channels>)
% A{105} = sprintf('%s',reinsert2);
% % Write cell A into txt
% fid = fopen(file, 'w');
% for i = 1:numel(A)
%     if A{i+1} == -1
%         fprintf(fid,'%s', A{i});
%         break
%     else
%         fprintf(fid,'%s\n', A{i});
%     end
% end
% mkdir('mean')
% cd('mean')
% movefile([folder3,'/',file],['meanspectrum_',num2str(c),'.spx']);
% clear ans fid i reinsert reinsert2 tline
% cd(folder3)
% 
% %%%%%%% PCAspectrum %%%%%%%
% reinsert=sprintf('%d,',pcaspectrum);
% reinsert2=['<Channels>',reinsert,'</Channels>'];
% % Read txt into cell A
% fid = fopen([folder,file],'r');
% i = 1;
% tline = fgetl(fid);
% A{i} = tline;
% while ischar(tline)
%     i = i+1;
%     tline = fgetl(fid);
%     A{i} = tline;
% end
% fclose(fid);
% % Change cell at position 105 (<Channels>)
% A{105} = sprintf('%s',reinsert2);
% % Write cell A into txt
% fid = fopen(file, 'w');
% for i = 1:numel(A)
%     if A{i+1} == -1
%         fprintf(fid,'%s', A{i});
%         break
%     else
%         fprintf(fid,'%s\n', A{i});
%     end
% end
% mkdir('PCA')
% cd('PCA')
% movefile([folder3,'/',file],['pcaspectrum_',num2str(c),'.spx']);
% clear ans fid i reinsert reinsert2 tline
% cd(folder3)
% 
% end
end
function [screen_int,facedata] = Cube_Generate(BinFile,isHex)
%MasterCube - Read a bin file from Dynamics & Build Interpolants

% This code is copyright Alex Foden and Ben Britton 09/04/2019
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
% OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
% NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
% BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN 
% ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
% CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

% Requirements:
% MATLAB R2018a or above
% MTEX version 5.2.beta2 or above
% Created by Alex Foden and Ben Britton 28/03/2019
% If you are using a CIF file not in the MTEX toolbox, you will need to add
% the full file path to the cif file to the phase file you are using

fileID = fopen(BinFile, 'r', 'ieee-le');
% if fileID == -1, error('Cannot open file: %s', filename); end
format = 'uint';
Data = fread(fileID, Inf, format);
fclose(fileID);
%find out the simulation resolution
cube_res=sqrt((size(Data,1)-12)/6)-1;

%build the face data for the cube
fd=zeros(cube_res+1,cube_res+1,6);
for n=1:6
    dstart=(n-1)*(cube_res+1)^2+1+2*n-1;
    dend=n*(cube_res+1)^2+2*n-1;
    fd(:,:,n)=flipud(rot90(reshape(Data(dstart:dend),cube_res+1,cube_res+1)));
end

%normalise
fd=fd-min(fd(:));
fd=fd/max(fd(:));
facedata=fd;
%sort this data
facedata(:,:,1)=fd(:,:,3); %x+
facedata(:,:,2)=fd(:,:,5); %y+
facedata(:,:,3)=fd(:,:,1); %z+

if isHex == 1 %add in a hexagonal fix
    facedata(:,:,4)=rot90(fd(:,:,1),2);
    facedata(:,:,5)=rot90(fd(:,:,2),2);
    facedata(:,:,6)=(fd(:,:,3));
else
    facedata(:,:,4)=fd(:,:,4); %x-
    facedata(:,:,5)=fd(:,:,6); %y-
    facedata(:,:,6)=fd(:,:,2); %z-
end

%build the interpolants
[gx,gy]=ndgrid(linspace(-1,1,cube_res+1),linspace(-1,1,cube_res+1));
screen_int.p1=griddedInterpolant(gx,gy,facedata(:,:,1),'cubic'); %x+
screen_int.p2=griddedInterpolant(gx,gy,facedata(:,:,2),'cubic'); %y+
screen_int.p3=griddedInterpolant(gx,gy,facedata(:,:,3),'cubic'); %z+
screen_int.p4=griddedInterpolant(gx,gy,facedata(:,:,4),'cubic'); %x-
screen_int.p5=griddedInterpolant(gx,gy,facedata(:,:,5),'cubic'); %y-
screen_int.p6=griddedInterpolant(gx,gy,facedata(:,:,6),'cubic'); %z-

screen_int.isHex=isHex;
  
% screen_int2.p1=griddedInterpolant(gx,gy,1+0*facedata(:,:,1),'cubic'); %x+
% screen_int2.p2=griddedInterpolant(gx,gy,2+0*facedata(:,:,2),'cubic'); %y+
% screen_int2.p3=griddedInterpolant(gx,gy,3+0*facedata(:,:,3),'cubic'); %z+
% screen_int2.p4=griddedInterpolant(gx,gy,4+0*facedata(:,:,4),'cubic'); %x-
% screen_int2.p5=griddedInterpolant(gx,gy,5+0*facedata(:,:,5),'cubic'); %y-
% screen_int2.p6=griddedInterpolant(gx,gy,6+0*facedata(:,:,6),'cubic'); %z-

end

function [i_data,screen_number] = Cube_Sample(xjs,yjs,zjs,screen_int,isHex)
%Sample a diffraction cube for intepolation

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

% inputs - 
% xjs, yjx, zjs = x y and z of sampling vectors
% screen_int = structure taken from MasterCube
% outputs -
% i_data = intensity
% screen_number = which screen was used from the cube

magv = sqrt(xjs.^2 + yjs.^2 + zjs.^2);
rjs=[xjs(:)./magv(:),yjs(:)./magv(:),zjs(:)./magv(:)];

if isHex == 1
    preA=[0.5*sqrt(2)/sqrt(3),1/sqrt(2),1/sqrt(3);...
        -sqrt(2)/sqrt(3),0,1/sqrt(3);...
        0.5*sqrt(2)/sqrt(3) -1/sqrt(2) 1/sqrt(3)];
    %rotate [111] onto [001] and put Y//a2
    rjs=rjs*preA';
end

%find which screen they should sample
[~,screen_number]=max(abs(rjs),[],2);
absmax = rjs(sub2ind(size(rjs),1:size(rjs,1),screen_number'));
signmax=sign(absmax);
screen_number(signmax==-1)=screen_number(signmax==-1)+3;

%sample the screens
i_data=zeros(size(screen_number));

%cycle through each screen
for num_screen=1:6
    s_inds=find(screen_number==num_screen);
    rjs_screen=rjs(s_inds,:)./transpose(repmat(absmax(s_inds).*signmax(s_inds),[3,1]));
    
    switch num_screen %these cases have been built through trialing
        case 1 %x+
            i_data(s_inds)=screen_int.p1(rjs_screen(:,2),-rjs_screen(:,3));
        case 2 %y+
             i_data(s_inds)=screen_int.p2(-rjs_screen(:,3),rjs_screen(:,1));
        case 3 %z+
            i_data(s_inds)=screen_int.p3(rjs_screen(:,2),rjs_screen(:,1));
        case 4 %x-
            i_data(s_inds)=screen_int.p4(rjs_screen(:,2),rjs_screen(:,3));
        case 5 %y-
            i_data(s_inds)=screen_int.p5(rjs_screen(:,3),rjs_screen(:,1));
        case 6 %z-
            i_data(s_inds)=screen_int.p6(-rjs_screen(:,2),rjs_screen(:,1));
    end
    

end

end


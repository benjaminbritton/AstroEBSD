function eang_out=conv_G_to_EA(g_in)
%CONV_G_TO_EA converts a gmatrix to Euler angles (rads)
%
%Follows the Euler Angle convention described in
%
% T.B. Britton, J. Jiang, Y. Guo, A. Vilalta-Clemente, D. Wallis,
% L. Hansen, A. Winkelmann, A.J. Wilkinson 
% Tutorial: Crystal orientations and EBSD — Or which way is up? 
% Materials Characterization (2016)
% http://dx.doi.org/10.1016/j.matchar.2016.04.008

%derived from
%eangs_in=[phi1 PHI phi2]
%R_x=@(theta)[1 0 0;0 cos(theta) sin(theta);0 -sin(theta) cos(theta)];
%R_z=@(theta)[cos(theta) sin(theta) 0;-sin(theta) cos(theta) 0;0 0 1];
%EBSD_Rot_in=R_z(eangs_in(3))*R_x(eangs_in(2))*R_z(eangs_in(1));

num_g=size(g_in,3);
eang_out=zeros(3,num_g);

for g=1:num_g
    gmat=g_in(:,:,g);
    
    if round(gmat(3,3),9) == 1
        eang=[atan2(gmat(1,2),gmat(1,1))/2,0,atan2(gmat(1,2),gmat(1,1))/2];
    else
        eang=[atan2(gmat(3,1),-gmat(3,2)),acos(gmat(3,3)),atan2(gmat(1,3),gmat(2,3))];
    end
    eang_out(:,g)=eang;
end
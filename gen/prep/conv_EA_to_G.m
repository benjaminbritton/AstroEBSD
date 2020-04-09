function gmat=conv_EA_to_G(eang)
%CONV_EA_TO_G - converts Euler angles (rads) to G matrix
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

% gmat=[   cos(eang(1))*cos(eang(3)) - cos(eang(2))*sin(eang(1))*sin(eang(3)), cos(eang(3))*sin(eang(1)) + cos(eang(2))*cos(eang(1))*sin(eang(3)), sin(eang(2))*sin(eang(3));
%              cos(eang(1))*sin(eang(3)) - cos(eang(2))*cos(eang(3))*sin(eang(1)), cos(eang(2))*cos(eang(1))*cos(eang(3)) - sin(eang(1))*sin(eang(3)), sin(eang(2))*cos(eang(3));
%                                                       sin(eang(2))*sin(eang(1)),                                -sin(eang(2))*cos(eang(1)),                      cos(eang(2))];

                                                  
G_eang1=[cos(eang(1)) sin(eang(1)) 0; -sin(eang(1)) cos(eang(1)) 0; 0 0 1];
G_eang2=[1 0 0; 0 cos(eang(2)) sin(eang(2)); 0 -sin(eang(2)) cos(eang(2))];
G_eang3=[cos(eang(3)) sin(eang(3)) 0; -sin(eang(3)) cos(eang(3)) 0; 0 0 1];
gmat=G_eang3*(G_eang2*G_eang1);
                                                  
                                                  
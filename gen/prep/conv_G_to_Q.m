function [ Q ] = conv_G_to_Q( G )
%conv_G_to_Q Converts G matrix into quaternions
%gives the quaternion as [w,XYZ]

G=G./det(G);

qtests=[0.5*sqrt(1 + G(1,1) + G(2,2) + G(3,3));
    0.5*sqrt(1 + G(1,1) - G(2,2) - G(3,3));
    0.5*sqrt(1 + G(2,2) - G(1,1) - G(3,3));
    0.5*sqrt(1 + G(3,3) - G(1,1) - G(2,2))];
[D,i]=max(qtests);
switch i
    case 1 %qw
        Qtemp(1)=D;
        Qtemp(2)=(G(3,2)-G(2,3))/(4*D);
        Qtemp(3)=(G(1,3)-G(3,1))/(4*D);
        Qtemp(4)=(G(2,1)-G(1,2))/(4*D);
    case 2 %qx
        Qtemp(1)=(G(3,2)-G(2,3))/(4*D);
        Qtemp(2)=D;
        Qtemp(3)=(G(1,2)+G(2,1))/(4*D);
        Qtemp(4)=(G(3,1)+G(1,3))/(4*D);
    case 3 %qw
        Qtemp(1)=(G(1,3)-G(3,1))/(4*D);
        Qtemp(2)=(G(1,2)+G(2,1))/(4*D);
        Qtemp(3)=D;
        Qtemp(4)=(G(2,3)+G(3,2))/(4*D);
    case 4 %qz
        Qtemp(1)=(G(2,1)-G(1,2))/(4*D);
        Qtemp(2)=(G(3,1)+G(1,3))/(4*D);
        Qtemp(3)=(G(2,3)+G(3,2))/(4*D);
        Qtemp(4)=D;
end

Qrot180=[cos(pi/2) sin(pi/2) 0 0]; % about x-axis, see Rollett lectures

Q=Q_MultArray(Qtemp',Qrot180');
function [ QOut_Array ] = Q_MultArray( Q1,Q2 )
%Q_MULT Multiply two quaternion arrays
%   TBB 2013 - VT edited

if size(Q1,2)==1 && size(Q2,2)==1
%     Q1 as (4,1) and Q2 as (4,1)
    QOut_Array=[Q1(1)*Q2(1) - Q1(2)*Q2(2) - Q1(3)*Q2(3) - Q1(4)*Q2(4);...
        Q1(1)*Q2(2) + Q1(2)*Q2(1) + Q1(3)*Q2(4) - Q1(4)*Q2(3);...
        Q1(1)*Q2(3) - Q1(2)*Q2(4) + Q1(3)*Q2(1) + Q1(4)*Q2(2);...
        Q1(1)*Q2(4) + Q1(2)*Q2(3) - Q1(3)*Q2(2) + Q1(4)*Q2(1)];
    
elseif size(Q1,2)== 1
%    Q1 as a (4,1) and Q2 (4,N)
QOut_Array=[Q1(1)*Q2(1,:) - Q1(2)*Q2(2,:) - Q1(3)*Q2(3,:) - Q1(4)*Q2(4,:);...
    Q1(1)*Q2(2,:) + Q1(2)*Q2(1,:) + Q1(3)*Q2(4,:) - Q1(4)*Q2(3,:);...
    Q1(1)*Q2(3,:) - Q1(2)*Q2(4,:) + Q1(3)*Q2(1,:) + Q1(4)*Q2(2,:);...
    Q1(1)*Q2(4,:) + Q1(2)*Q2(3,:) - Q1(3)*Q2(2,:) + Q1(4)*Q2(1,:)];
  
elseif size(Q1,2)> 1 && size(Q1,1)== 4 
%  Q1 as a (4,N) and Q2 (4,N)
  QOut_Array=[Q1(1,:).*Q2(1,:) - Q2(2,:).*Q1(2,:) - Q2(3,:).*Q1(3,:) - Q2(4,:).*Q1(4,:);...
        Q1(1,:).*Q2(2,:) + Q1(2,:).*Q2(1,:) + Q1(3,:).*Q2(4,:) - Q1(4,:).*Q2(3,:);...
        Q1(1,:).*Q2(3,:) + Q2(1,:).*Q1(3,:) - Q1(2,:).*Q2(4,:) + Q2(2,:).*Q1(4,:);...
        Q1(1,:).*Q2(4,:) + Q1(2,:).*Q2(3,:) + Q2(1,:).*Q1(4,:) - Q1(3,:).*Q2(2,:)];
    
else
    
    %should extend to a symmetric version with an extra end index
    %i.e. (Y,X,[Qw,Qx,Qy,Qz],symmetry)
    
%   Array must be (Y,X,[Qw,Qx,Qy,Qz])
% QOut_Array_w=Q1(:,:,1).*Q2(:,:,1) - Q1(:,:,2).*Q2(:,:,2) - Q1(:,:,3).*Q2(:,:,3) - Q1(:,:,4).*Q2(:,:,4);
% QOut_Array_x=Q1(:,:,1).*Q2(:,:,2) + Q1(:,:,2).*Q2(:,:,1) + Q1(:,:,3).*Q2(:,:,4) - Q1(:,:,4).*Q2(:,:,3);
% QOut_Array_y=Q1(:,:,1).*Q2(:,:,3) - Q1(:,:,2).*Q2(:,:,4) + Q1(:,:,3).*Q2(:,:,1) + Q1(:,:,4).*Q2(:,:,2);
% QOut_Array_z=Q1(:,:,1).*Q2(:,:,4) + Q1(:,:,2).*Q2(:,:,3) - Q1(:,:,3).*Q2(:,:,2) + Q1(:,:,4).*Q2(:,:,1);

P1_Q=repmat(Q1(:,:,1),[1,1,4]).*Q2;
P2_Q=repmat(Q1(:,:,2),[1,1,4]).*Q2;
P3_Q=repmat(Q1(:,:,3),[1,1,4]).*Q2;
P4_Q=repmat(Q1(:,:,4),[1,1,4]).*Q2;

QOut_Array_w=P1_Q(:,:,1)-P2_Q(:,:,2)-P3_Q(:,:,3)-P4_Q(:,:,4);
QOut_Array_x=P1_Q(:,:,2)+P2_Q(:,:,1)+P3_Q(:,:,4)-P4_Q(:,:,3);
QOut_Array_y=P1_Q(:,:,3)-P2_Q(:,:,4)+P3_Q(:,:,1)+P4_Q(:,:,2);
QOut_Array_z=P1_Q(:,:,4)+P2_Q(:,:,3)-P3_Q(:,:,2)+P4_Q(:,:,1);

% QOut_Array_Norm=sqrt(QOut_Array_w.^2+QOut_Array_x.^2+QOut_Array_y.^2+QOut_Array_z.^2);

QOut_Array=zeros(size(Q1));
QOut_Array(:,:,1)=QOut_Array_w;%./QOut_Array_Norm;
QOut_Array(:,:,2)=QOut_Array_x;%./QOut_Array_Norm;
QOut_Array(:,:,3)=QOut_Array_y;%./QOut_Array_Norm;
QOut_Array(:,:,4)=QOut_Array_z;%./QOut_Array_Norm;

end
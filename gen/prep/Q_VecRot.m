function [ output_args ] = Q_VecRot( quat,vec )
%QUAT_VEC_ROT Summary of this function goes here
%   Detailed explanation goes here

q2_v2=quat(2,:).*vec(2);
q3_v1=quat(3,:).*vec(1);
q3_v3=quat(3,:).*vec(3);
q4_v2=quat(4,:).*vec(2);
q2_v3=quat(2,:).*vec(3);
q4_v1=quat(4,:).*vec(1);

output_args=([ vec(1) + 2.*quat(3,:).*(q2_v2 - q3_v1) + 2.*quat(1,:).*(q3_v3 - q4_v2) + 2.*quat(4,:).*(q2_v3 - q4_v1);
    vec(2) - 2.*quat(2,:).*(q2_v2 - q3_v1) - 2.*quat(1,:).*(q2_v3 - q4_v1) + 2.*quat(4,:).*(q3_v3 - q4_v2); 
    vec(3) + 2.*quat(1,:).*(q2_v2 - q3_v1) - 2.*quat(2,:).*(q2_v3 - q4_v1) - 2.*quat(3,:).*(q3_v3 - q4_v2)]);


end


function [gmat] = secular_equation(kmat,s1,s2,s1s,s2s,gamma,gammas,h)
% Description:
%   This function calculates the determinant of the characteristic matrix
%
% Outputs:
%   gmat - determinants
%
% Author:
%   Yuxuan Jiang, @ Wellman Center for Photomedicine
%   2025-12
%
% Github: https://github.com/JanYxuan/Bilayer-wave-dispersion.git
%
% Reference:
%   [1] X Feng, GY Li, Y Jiang, O Shortt-Nguyen, SH Yun. Optical coherence
%   elastography measures mechanical tension in the lens and capsule. Acta
%   Biomaterialia, 199:252-261, 2025.

%% ------ Components ------
seq=ones(1,length(kmat));
%
m11=s1;
m12=s2;
m13=-s1;
m14=-s2;
m15=-s1s;
m16=-s2s;
%
m21=1*seq;
m22=1*seq;
m23=1*seq;
m24=1*seq;
m25=-1*seq;
m26=-1*seq;
%
m31=gamma*(1+s1.^2);
m32=gamma*(1+s2.^2);
m33=gamma*(1+s1.^2);
m34=gamma*(1+s2.^2);
m35=-gammas*(1+s1s.^2);
m36=-gammas*(1+s2s.^2);
%
m41=gamma*s1.*(1+s2.^2);
m42=gamma*s2.*(1+s1.^2);
m43=-gamma*s1.*(1+s2.^2);
m44=-gamma*s2.*(1+s1.^2);
m45=-gammas*s1s.*(1+s2s.^2);
m46=-gammas*s2s.*(1+s1s.^2);
%
m51=(1+s1.^2).*exp(s1.*kmat*h);
m52=(1+s2.^2).*exp(s2.*kmat*h);
m53=(1+s1.^2).*exp(-s1.*kmat*h);
m54=(1+s2.^2).*exp(-s2.*kmat*h);
m55=0;
m56=0;
%
m61=s1.*(1+s2.^2).*exp(s1.*kmat*h);
m62=s2.*(1+s1.^2).*exp(s2.*kmat*h);
m63=-s1.*(1+s2.^2).*exp(-s1.*kmat*h);
m64=-s2.*(1+s1.^2).*exp(-s2.*kmat*h);
m65=0;
m66=0;

%% --------- Characteristics matrix ----------
M=zeros(6,6,length(kmat));
M(1,1,:)=m11; M(1,2,:)=m12; M(1,3,:)=m13; M(1,4,:)=m14; M(1,5,:)=m15; M(1,6,:)=m16;
M(2,1,:)=m21; M(2,2,:)=m22; M(2,3,:)=m23; M(2,4,:)=m24; M(2,5,:)=m25; M(2,6,:)=m26; 
M(3,1,:)=m31; M(3,2,:)=m32; M(3,3,:)=m33; M(3,4,:)=m34; M(3,5,:)=m35; M(3,6,:)=m36;
M(4,1,:)=m41; M(4,2,:)=m42; M(4,3,:)=m43; M(4,4,:)=m44; M(4,5,:)=m45; M(4,6,:)=m46;
M(5,1,:)=m51; M(5,2,:)=m52; M(5,3,:)=m53; M(5,4,:)=m54; M(5,5,:)=m55; M(5,6,:)=m56;
M(6,1,:)=m61; M(6,2,:)=m62; M(6,3,:)=m63; M(6,4,:)=m64; M(6,5,:)=m65; M(6,6,:)=m66;

% ----- determinant ------
gmat=zeros(1,length(kmat));
for ii=1:length(gmat)
    Mt=M(:,:,ii);
    g=det(Mt);
    gmat(ii)=g;
end


end


function [ imin,kmat,gmat3,gmat2,flag_mark ] = SearchRoots( w,cmat,kim,alpha,beta,gamma,alpha_sub,beta_sub,gamma_sub,rho,rhosub,h,trigger )
% Description:
%   This function find the local minimum of the determinants
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

kre=w./cmat;
[KRE,KIM]=meshgrid(kre,kim);
kmat=KRE-1i*KIM;
kmat2=reshape(kmat,1,[]);
w2kmat=w./kmat2;

delta=(2*beta-rho*w2kmat.^2).^2-4*gamma*(alpha-rho*w2kmat.^2);
s1=sqrt(((2*beta-rho*w2kmat.^2)+sqrt(delta))/(2*gamma));
s2=sqrt(((2*beta-rho*w2kmat.^2)-sqrt(delta))/(2*gamma));
delta_sub=(2*beta_sub-rhosub*w2kmat.^2).^2-4*gamma_sub*(alpha_sub-rhosub*w2kmat.^2);
s1_sub=sqrt(((2*beta_sub-rhosub*w2kmat.^2)+sqrt(delta_sub))/(2*gamma_sub));
s2_sub=sqrt(((2*beta_sub-rhosub*w2kmat.^2)-sqrt(delta_sub))/(2*gamma_sub));
if trigger==1
    s2_sub=-s2_sub;
end

gmat=secular_equation(kmat2,s1,s2,s1_sub,s2_sub,gamma,gamma_sub,h);
gmat2=log10(abs(gmat));
gmat3=reshape(gmat2,length(kim),length(kre));
[~,~,~,imin]=extrema2(gmat3);
[imin,flag_mark]=eliminateData(gmat3,imin);

end

%% ========= FUNCTION =========
function [ind,ind1,ind2,ind3,ind4]=findMargin(gmat2)
[m,n]=size(gmat2);
ind1=2:m; % left,kre=0
ind2=m*(n-1)+2:m*n; % right, kre=max
ind3=m:m:m*n; % up
ind4=1:m:m*(n-1)+1; %down
ind=[ind1,ind2,ind3,0];
ind=unique(ind);
end

function [imin2,flag_mark]=eliminateData(gmat2,imin)
[ind,ind1,ind2,ind3,ind4]=findMargin(gmat2);
a1=intersect(ind,imin);
imin2=setdiff(imin,a1);
flag_mark=DetermineEdge(imin,ind1,ind2,ind3,ind4);
end

function [flag_mark]=DetermineEdge(imin,ind1,ind2,ind3,ind4)
if isempty(imin)
    flag_mark=-1; % nan
elseif length(imin)>=2
    flag_mark=5; % more than two roots
elseif length(imin)==1&&ismember(imin,ind1)
    flag_mark=1; % left
elseif length(imin)==1&&ismember(imin,ind2)
    flag_mark=2; % right
elseif length(imin)==1&&ismember(imin,ind3)
    flag_mark=3; % up
elseif length(imin)==1&&ismember(imin,ind4)
    flag_mark=4; % down
else
    flag_mark=0; % center
end
end
function [phase_vel] = SolveStressedBilayerWaveDispersion(paras,freq_range)
% Description:
%  This function adaptively determines the search range in the k-domain, and 
%  identifies local minima within the k-domain.
%
% Inputs:
%   paras - material and geometric parameters of the bilayer system
%   freq_range - the frequency points specified by the user
%
% Outputs:
%   phase_vel - phase velocities at the specified frequency points
%
% Author:
%   Yuxuan Jiang, @ Wellman Center for Photomedicine
%   2025-12
%   jiangyx96@gmail.com
%
% Github: https://github.com/JanYxuan/Bilayer-wave-dispersion.git
%
% Reference:
%   [1] X Feng, GY Li, Y Jiang, O Shortt-Nguyen, SH Yun. Optical coherence
%   elastography measures mechanical tension in the lens and capsule. Acta
%   Biomaterialia, 199:252-261, 2025.

[mu,sigma,mu_s,sigma2,rho,rhosub,h]=deal(paras(1),paras(2),paras(3),paras(4),paras(5),paras(6),paras(7));
fmat = UpdataFreqList(freq_range);

% --- Initialization ---
ct2=sqrt((mu_s+sigma2/3)/rhosub);
dc=2e-3;
n1=length(fmat);
wmat=2*pi*fmat;

% --- inc paras ---
alpha=mu+sigma/3;
gamma=mu-sigma*2/3;
beta=mu-sigma/6;
if gamma<=0
    error('\mu < 1.5\sigma');
end
alpha_sub=mu_s+sigma2/3;
gamma_sub=mu_s-sigma2*2/3;
beta_sub=mu_s-sigma2/6;
if gamma_sub<=0
    error('\mu_s < 1.5\sigma_s');
end

% -- Settings --
ctag_mat=zeros(1,length(wmat));
kim_mat=zeros(1,length(wmat));
cgap0=0.5;
kimgap0=40;
flag=4;
trigger=0;

% -- Scan each frequency --
for ii=1:n1
    w=wmat(ii);
    if ii==1
        c01=ct2;
        kim01=0;
    else
        c01=ctag_mat(ii-1);
        kim01=kim_mat(ii-1);
    end
    if ii<=2
        cgap=cgap0;
        kimgap=kimgap0;
    else
        cgap1=abs(ctag_mat(ii-1)-ctag_mat(ii-2));
        kimgap1=abs(kim_mat(ii-1)-kim_mat(ii-2));
        if cgap1<dc
            cgap1=dc;
        end
        if kimgap1<1e-2
            kimgap1=kimgap0;
        end
        cgap=min(cgap1*3,cgap0);
        kimgap=max(kimgap1*3,kimgap0);
    end
    if ii<=10 && all(kim_mat==0)
        cmat=0.8*ct2:dc:1.3*ct2;
    else
        clow=max(dc,c01-cgap*1);
        cmat=clow:dc:c01+cgap*3;
    end
    if ii<=10
        kim=linspace(0,max(w./cmat)/5,50);
    else
        kimlow=max(0,kim01-kimgap*2);
        kim=linspace(kimlow,kim01+kimgap*4,50);
    end
    
    [imin,kmat,~,gmat2,flag]=SearchRoots(w,cmat,kim,alpha,beta,gamma,alpha_sub,beta_sub,gamma_sub,rho,rhosub,h,trigger);
    if flag==5&&isempty(imin)
        flag=-1;
    end

    % ---- Resampling for the failed case ---
    if flag==1
        clow=max(dc,c01-cgap*10);
        cmat2=clow:dc:c01+cgap/2;
        [imin,kmat,~,gmat2,flag]=SearchRoots(w,cmat2,kim,alpha,beta,gamma,alpha_sub,beta_sub,gamma_sub,rho,rhosub,h,trigger);
        
    elseif flag==2
        clow=max(dc,c01-cgap*0.1);
        cmax=max(1.2*ct2,c01+cgap*5);
        cmat2=clow:dc:cmax;
        [imin,kmat,~,gmat2,flag]=SearchRoots( w,cmat2,kim,alpha,beta,gamma,alpha_sub,beta_sub,gamma_sub,rho,rhosub,h,trigger);
        
    elseif flag==3
        kimlow=max(0,kim01-kimgap*1);
        kimhigh=kim01+max(kimgap1*5,2000);
        kim2=linspace(kimlow,kimhigh,50);
        [imin,kmat,~,gmat2,flag]=SearchRoots( w,cmat,kim2,alpha,beta,gamma,alpha_sub,beta_sub,gamma_sub,rho,rhosub,h,trigger);
    
    elseif flag==-1
        clow=max(dc,c01-cgap*10);
        cmat2=clow:dc:c01+cgap*20;
        [imin,kmat,~,gmat2,flag]=SearchRoots(w,cmat2,kim,alpha,beta,gamma,alpha_sub,beta_sub,gamma_sub,rho,rhosub,h,trigger);
    end
    
    if flag==5
        % length of (imin)>=2
        if ii==1
            ctmp=w./real(kmat(imin)); pos=find(ctmp<ct2);
            if length(pos)>1
                tmp=gmat2(imin); [~,pos]=min(tmp);
            end
            imin=imin(pos);
        elseif ii>=2&&ii<=3
            tmp=gmat2(imin); [~,pos]=min(tmp); imin=imin(pos);
        else
            ctmp=w./real(kmat(imin));
            cesti=2*ctag_mat(ii-1)-ctag_mat(ii-2);
            [~,pos]=min(abs(ctmp-cesti)); imin=imin(pos);
        end
    end
    
    ktag=kmat(imin);
    ktag_re=real(ktag);
    ktag_im=imag(ktag);
    ctag=w/ktag_re;
    ctag_mat(ii)=ctag;
    kim_mat(ii)=abs(ktag_im);
    
    if any(ctag_mat>=ct2*1.02)
        trigger=1;
    else
        trigger=0;
    end
end

phase_vel=interp1(fmat,ctag_mat,freq_range,'linear','extrap');

end


function [f_out] = UpdataFreqList(f_user)
default_step = 200;
fmin = min(f_user);
fmax = max(f_user);
if fmin < 1
    f_start = fmin;
else
    f_start = 1;
end

if fmax < 1
    error('max freq should be no less than 1 Hz');
end
if fmax < default_step
    n_points = max(10, ceil((fmax - f_start) + 1));
    f_out = linspace(f_start, fmax, n_points);
else
    f_out = f_start:default_step:fmax;
    if f_out(end) < fmax
        f_out(end+1) = f_out(end) + default_step;
    end
end
f_out = unique(f_out);
end

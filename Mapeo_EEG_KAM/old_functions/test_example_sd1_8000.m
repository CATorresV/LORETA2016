close all
clear 
clc
%Test 2
% SNRs=[14 12 7 0 -5];
SNRs = 7;
kS=length(SNRs);%number of SNR
nTRIAL=30; %number of trials
Nas=3; %number of active sources
methods={'LOR','SPM','IRA-TI-L1','IRA-TI-L2'};
Nmeth=length(methods);%number of methods
Residual_error=zeros(nTRIAL,kS,Nmeth);
Standard_error=zeros(nTRIAL,kS,Nmeth);
TAIy=zeros(nTRIAL,kS,Nmeth);
SAIy=zeros(nTRIAL,kS,Nmeth);
EMDy=zeros(nTRIAL,kS,Nmeth);
Nd=8000;
switch Nd
    case 20000
        load Model_head
        vert = vol.vert;
        face = vol.face;
        load L_neigbor
        load M_1020
    case 8000
        load L
        M = L; clear L;
        load QG
        load head_model
end
nf = [6000 2200 4300]%randi([1 size(vert,1)],1,Nas);
%nf=5534
RealX=zeros(length(QG),nTRIAL,kS);
EstX=zeros(length(QG),nTRIAL,kS,Nmeth);
%Dynamic Constrained Inverse Problem
%%
%EMD
D = full(QG);
D = D-min(min(D));
D = D./max(max(D));
D = 1-D;
extra_mass_penalty= -1;
%%
%Spatial basis
[M_s,Phi_s]=spatial_basis(Nd);
[d,s]=size(M_s);
%s:number of states
%d:number of measurements%%
%Reduction of computational load
[U,S,V] = svds(M_s,d);
d_r=d;
% s_cum=cumsum(diag(S))/sum(diag(S));
% ind_s=find(s_cum>0.999);
% d_r=ind_s(1)
% [U,S,V] = svds(M_s,d_r);
%SNR (dB)
%load DataTest
kt=1
for SNR=SNRs,
    SNR%%
    for TRIAL=1:nTRIAL,
        %%
        TRIAL
        %Total Samples
        N=250;
        %Simulation time (seconds)
        inter=[0,1];
        %SNR (dB)
        %EEG Simulation
        [y,x,t]=eeg_dynamic_simulation_snr(N,inter,SNR,Nd,Nas,nf);
        %y=y./max(max(abs(y)));
        RealX(:,TRIAL,kt)=mean(x.^2,2);
%         Plot3DBrain(RealX(:,TRIAL),vert,face,QG);
%         pause
        [K]=length(t);
        %K:total number of samples
        for meth=1:Nmeth,
            methods{meth}
            switch meth,
                case 1,
                    %%
                    % --------------------------
                    inverse.type = 'LOR'; %Loreta
                    % inverse.type = 'EBB'; %Beamformer
                    % inverse.type = 'GS'; %MSP
                    J = spm_eeg_invert_red(y,M,QG,vert,face,inverse);
                    x_est=J;
                    y_est=M*x_est;
                    %%
                case 2,
                    %%
                    %inverse.type = 'LOR'; %Loreta
                    % inverse.type = 'EBB'; %Beamformer
                    inverse.type = 'GS'; %MSP
                    
                    J = spm_eeg_invert_red(y,M,QG,vert,face,inverse);
                    x_est=J;
                    y_est=M*x_est;
           
                case 3,
                    [UU,ss,vv]    = spm_svd((M*M'),exp(-16));
                    AA     = UU';                 % spatial projector A
                    UM    = AA*M;
                    Mold = M;
                    M = UM;
                    yold = y;
                    y = AA*y;
                    M_s_old = M_s;
                    M_s = M*Phi_s;
                    
                    c_est=zeros(s,K);
                    rp_est=zeros(2,K);
                    %a priori estimate
                    c1=ones(s,1);
                    rp1=[1.5 0.5];
                    %Iterative Regularized Solution
                    for k=1:K
                        k
                        %Parameter selection by Generalized Cross Validation
                         if k<=1
                            En = sum(y.^2,1);
                            [~,ind] = max(En);
                            rp=GCV_multi_IRA_TI_L1(M_s,y(:,ind),rp1,c1)
                        else
                            rp=rp1;
                        end
                        rp_est(:,k)=rp(:);
                        rho=rp(1);
                        lambda=rp(2);
                        %Iterative Regularized Algorithm for Inverse Problem Solution
                        %State estimation
                        W=diag(abs(c1));
                        try   
                            c_est(:,k)= pinv(W*M_s'*M_s+lambda^2*W+rho^2*eye(s))*W*(M_s'*y(:,k)+lambda^2*c1);
                        catch
                            c_est(:,k)= pinvecon(W*M_s'*M_s+lambda^2*W+rho^2*eye(s))*W*(M_s'*y(:,k)+lambda^2*c1);
                        end
                        %store past iterations
                        c1=c_est(:,k);
                        rp1=rp;
                    end
                    x_est=Phi_s*c_est;
                    y_est=AA'*M_s*c_est;
                    
                    y = yold;
                    M_s = M_s_old;
                    M = Mold;
                case 4,
                     [UU,ss,vv]    = spm_svd((M*M'),exp(-16));
                    AA     = UU';                 % spatial projector A
                    UM    = AA*M;
                    Mold = M;
                    M = UM;
                    yold = y;
                    y = AA*y;
                    M_s_old = M_s;
                    M_s = M*Phi_s;
                    [U,S,V] = svds(M_s,d);

                    c_est=zeros(s,K);
                    rp_est=zeros(2,K);
                    %a priori estimate
                    c1=ones(s,1);
                    rp1=[0.1 0.1];
                    %Iterative Regularized Solution
                    for k=1:K
                        k
                        %Parameter selection by Generalized Cross Validation
                         if k<=1
                            En = sum(y.^2,1);
                            [~,ind] = max(En);
                            rp=GCV_multi_IRA_TI_L2_svd(M_s,U,S,V, y(:,ind), rp1,c1)
                        else
                            rp=rp1;
                        end
                        rp_est(:,k)=rp(:);
                        rho=rp(1);
                        lambda=rp(2);
                        %Iterative Regularized Algorithm for Inverse Problem Solution
                        %State estimation
                       c_est(:,k)= V*diag(1./(diag(S).^2+(lambda^2+rho^2)))*V'*(M_s'*y(:,k)+lambda^2*c1);
                       %store past iterations
                        c1=c_est(:,k);
                        rp1=rp;
                    end
%                     RP
                    x_est=Phi_s*c_est;
                    y_est=AA'*M_s*c_est;
                    
                    y = yold;
                    M_s = M_s_old;
                    M = Mold;
                    
            end
            
            %%
            Residual_error(TRIAL,kt,meth)=norm(y-y_est,'fro')^2;
            Standard_error(TRIAL,kt,meth)=norm(x-x_est,'fro')^2;
            SSerr=Residual_error(TRIAL,kt,meth);
            SStot=norm(detrend(y),'fro')^2;
            TAIy(TRIAL,kt,meth)=1-SSerr/SStot;
            SAIy(TRIAL,kt,meth)=error_sai(x,x_est, QG,vert,face);
            EMDy(TRIAL,kt,meth)=emd_hat_gd_metric_mex(mean(x.^2,2),mean(x_est.^2,2),D,extra_mass_penalty)
            EstX(:,TRIAL,kt,meth)=mean(x_est.^2,2);
        end
    end
    kt=kt+1;
end
%%
save results_example_sd1_8000 RealX Residual_error Standard_error TAIy SAIy EstX EMDy
%%
TRIAL=1
Plot3DBrain(RealX(:,TRIAL),vert,face,QG);
Plot3DBrain(EstX(:,TRIAL,1),vert,face,QG);
Plot3DBrain(EstX(:,TRIAL,2),vert,face,QG);
Plot3DBrain(EstX(:,TRIAL,3),vert,face,QG);
Plot3DBrain(EstX(:,TRIAL,4),vert,face,QG);

figure
plot(t,y+diag(cumsum(max(abs(y(1:end,:)'))))*ones(size(y)))
h1=xlabel('Time(s)');
set(h1,'Interpreter','latex','fontsize',16);
h1=ylabel('Amplitude(V)');
set(gca,'YTick',[]);
set(h1,'Interpreter','latex','fontsize',16);
h1=title('EEG');
set(h1,'Interpreter','latex','fontsize',16);


%%
for i=1:size(EMDy,3)
    aux{i} = squeeze(EMDy(:,:,i));
end

figure('Position',[300 300 500 250])
h = [1 0 0; 0 1 0; 0 0 1];
aboxplot(aux,'labels',SNRs);
xlabel('SNR (dB)');
ylabel('EMD');
legend('Lor','MSP','IRA-L1','IRA-L2')
%%
for i=1:size(TAIy,3)
    aux{i} = squeeze(TAIy(:,:,i));
end

figure('Position',[300 300 500 250])
h = [1 0 0; 0 1 0; 0 0 1];
aboxplot(aux,'labels',SNRs);
xlabel('SNR (dB)');
ylabel('TAIy');
legend('Lor','MSP','IRA-L1','IRA-L2')


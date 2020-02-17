%Grupo de automática
%Mapping testing
%Andres Marino Alvarez Meza 
%Cristian Alejandro Torres 2016
%% function paths
clear all; close all;clc
addpath(genpath('Mapeo_EEG_KAM'))
addpath(genpath('Aux_Mat_AM_Mapping'))
%% Read transformation matrix model (standard):
load L %lead field matrix -> conductivity matrix
M = full(L); clear L
load QG %spatial relationships
load head_model %brain model isotropic
%% Parameters for simulating EEG
inter = [0 1];%time interval
Fs = 250; %sample frequency
t = linspace(inter(1),inter(2),Fs);
%% SPM configurations -> method - decomposition - Covariance matrix
nreps=[100 2500 3500 5000 8000];%round(linspace(100,8000,10));%
Ndip = 1;%# of simulated sources
dipamp = [1500; 2000; 1000; 500].*1e-9; %dip amplitudes
dipfreq = [20 15 10 25];%dip frequencies
SNR = [ 3 7 10 ];%dB
%methods
%{['MMN','L','L'],['MMN','L','G'],['MMN','G','L'],['MMN','G','G'],...
method = {['LOR','L','L'],['LOR','L','G'],['LOR','G','L'],['LOR','G','G']... %,...%['MSP','L','L'],['MSP','L','G'],['MSP','G','L'],['MSP','G','G'],...
          ['EBB','L','L'],['EBB','L','G'],['EBB','G','L'],['EBB','G','G']}%,...
          %['SFL','L','L'],['SFM','L','L'],['SFB','L','L']};
Error = zeros(length(nreps),length(method));
%% main loop
c = 1;
for i=1:length(nreps)
    for j = 1 : length(SNR)
        close all;
        fprintf('Simulating EEG...')
        [y,x] = eeg_stationary_simulationA(t,SNR(j),Ndip,M,QG,vert,face,dipamp,dipfreq,nreps(i)); %100,5000,8000
        fprintf('done\n')
        e0 = sum(x.^2,2);
        e0 = e0/max(e0);
        %% mapping method
        for z = 1 : length(method)
            tic;
            fprintf('It:%d/%d',c,nreps*length(SNR)*length(method));
            inverse.type = method{z}(1:3); %inverse method
            inverse.type_cov = method{z}(4); %covariancefunction
            inverse.type_eig = method{z}(5); %decomposition type
            %mapping
            if strcmp(inverse.type,'SFL')
                inverse.type='LOR';
                [J,AYYA] = spm_eeg_invert_sinF(y,M,QG,vert,face,inverse);
                em = sum(J.^2,2);
                em = em/max(em);
                [e1{j}(i,z),e2{j}(i,z),e3{j}(i,z)] = DLE(e0, em, vert, Ndip);
            elseif strcmp(inverse.type,'SFM')
                inverse.type='MSP';
                [J,AYYA] = spm_eeg_invert_sinF(y,M,QG,vert,face,inverse);
                em = sum(J.^2,2);
                em = em/max(em);
                [e1{j}(i,z),e2{j}(i,z),e3{j}(i,z)] = DLE(e0, em, vert, Ndip);
            elseif strcmp(inverse.type,'SFB')
                inverse.type='EBB';
                [J,AYYA] = spm_eeg_invert_sinF(y,M,QG,vert,face,inverse);
                em = sum(J.^2,2);
                em = em/max(em);
                [e1{j}(i,z),e2{j}(i,z),e3{j}(i,z)] = DLE(e0, em, vert, Ndip);
            else
                [J,AYYA] = spm_eeg_invert_AM(y,M,QG,vert,face,inverse);
                em = sum(J.^2,2);
                em = em/max(em);
                [e1{j}(i,z),e2{j}(i,z),e3{j}(i,z)] = DLE(e0, em, vert, Ndip);
            end
            
             c= c+1;
             tc(c) = toc;
             fprintf('etime %.2f\n',tc(i))
        end
    end
end
save('Results_20_06_2016_conc_10','e1','e2','e3')
%% plot boxplots
close all;clc
for j = 1 : length(SNR)
    figure
    boxplot(e1{j})
    set(gca,'XTickLabel',method)
     xticklabel_rotate(1:length(method),45,method)
    title(['e1 SNR=' num2str(SNR(j))])
    
    figure
    boxplot(e2{j})
     set(gca,'XTickLabel',method)
    xticklabel_rotate(1:length(method),45,method)
    title(['e2 SNR=' num2str(SNR(j))])
     
    figure
    boxplot(e3{j})
     set(gca,'XTickLabel',method)
    xticklabel_rotate(1:length(method),45,method)
    title(['e3 SNR=' num2str(SNR(j))])
     
end
showfigs_c(4)

%%

        
        
%         %% Ploting signals
% %         fprintf('Plotting EEG...')
% %         figure
% %         subplot(4,2,1)
% %         plot(t,x,'b')
% %         legend({'Clean EEG'})
% %         subplot(4,2,3)
% %         plot(t,y,'r')
% %         legend({['Noisy EEG SNR=' num2str(SNR) '[dB]']})
% %         subplot(4,2,2)
% %         plot(t,y_n,'r')
% %         legend({'My Loretta solution'})
% %         subplot(4,2,5)
% %         plot(t,Y1,'r')
% %         legend({'Minimum Norm'})
% %         
% %         subplot(4,2,4)
% %         plot(t,Y0,'r')
% %         legend({'Loretta solution'})
% %         subplot(4,2,6)
% %         plot(t,Y2,'r')
% %         legend({'EBB solution'})
% %         subplot(4,2,7)
% %         plot(t,Yk,'r')
% %         legend({'Kernel based'})
%         %% Source plotting
%         
%         D = full(QG);
%         D = D-min(min(D));
%         D = D./max(max(D));
%         D = 1-D;
%         extra_mass_penalty= -1;
%         
% %         figure
% %         subplot(611)
%         ex = sum(x.^2,2);
%         ex = ex/max(ex);
% %         stem(ex);
% %         title('Original')
%         
%         
% %         subplot(612)
%         e0 = sum(X_s.^2,2);
%         e0 = e0/max(e0);
% %         stem(e0);
% %         title(['my LOR'])
%         [er0_1,er0_2,ep0]=DLE(ex,e0,vert,Ndip);
%         
% %         subplot(613)
%         e1 = sum(J0.^2,2);
%         e1 = e1/max(e1);
% %         stem(e1);
% %         title(['LOR - ' inverse1.type_cov])
%         [er1_1,er1_2,ep1]=DLE(ex,e1,vert,Ndip);
%         
%         
%         
% %         subplot(614)
%         e2 = sum(J1.^2,2);
%         e2 = e2/max(e2);
% %         stem(e2);
% %         title(['MMN - ' inverse1.type_cov])
%         [er2_1,er2_2,ep2]=DLE(ex,e2,vert,Ndip);
%         
% %         subplot(615)
%         e3 = sum(J2.^2,2);
%         e3 = e3/max(e3);
% %         stem(e3);
% %         title(['EBB - ' inverse1.type_cov])
%         [er3_1,er3_2,ep3]=DLE(ex,e3,vert,Ndip);
%         
% %         subplot(616)
%         ek = sum(Jk.^2,2);
%         ek = ek/max(ek);
% %         stem(ek);
% %         title(['kernel - ' inverse1.type_cov])
%         [erk_1,erk_2,ep4]=DLE(ex,ek,vert,Ndip);
%         
%         if graf3==1
%             %plot source energies on the brain model
%             Plot3DBrain(ex,vert,face,QG); colorbar
%             title('Clean source')
%             
%             Plot3DBrain(e0,vert,face,QG); colorbar
%             reo = norm(ex-e0);
%             title(['my LOR  ' 'MSE=' num2str(reo)])
%             
%             Plot3DBrain(e1,vert,face,QG); colorbar
%             reo = norm(ex-e1);
%             title(['LOR - ' inverse1.type_cov 'MSE=' num2str(reo)])
%             
%             Plot3DBrain(e2,vert,face,QG); colorbar
%             reo = norm(ex-e2);
%             title(['MMN - ' inverse1.type_cov 'MSE=' num2str(reo)])
%             
%             Plot3DBrain(e3,vert,face,QG); colorbar
%             reo = norm(ex-e3);
%             title(['EBB - ' inverse1.type_cov 'MSE=' num2str(reo)])
%             
%             
%             Plot3DBrain(ek,vert,face,QG); colorbar
%             reo = norm(ex-ek);
%             title(['Kernel - ' inverse1.type_cov 'MSE=' num2str(reo)])
%         else
%             reo0 = norm(ex-e0);
%             reo1 = norm(ex-e1);
%             reo2 = norm(ex-e2);
%             reo3 = norm(ex-e3);
%             reok = norm(ex-ek);
%             cade = 'e0 = Lor, e1 = LOr AM, e2 = MMN, e3 = EBB, ek = kernel'; 
%         end
%         er_cad = 'euclidean, norml1, percent';
%         Error(i,:) = [er0_1 er1_1 er2_1 er3_1 erk_1 er0_2 er1_2 er2_2 er3_2 erk_2 ep0 ep1 ep2 ep3 ep4 reo0 reo1 reo2 reo3 reok];
%         %save([ruta day 'solutionEEG_recons_Ndip' num2str(Ndip) '_SNR' num2str(SNR) '_it' num2str(i) '.mat'],'x','y','cad','X_s','J0','J1','J2','Jk');
%         
%         
%     end
%     %save([ruta day 'error_50_SNR' num2str(SNR) '_Ndip' num2str(Ndip) '.mat'],'er_cad','Error');
% end
% 
% %%
% 
% figure
% boxplot(Error)

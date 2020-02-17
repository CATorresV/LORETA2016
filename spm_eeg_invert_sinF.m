function [J,AYYA] = spm_eeg_invert_sinF(Y,L,QG,vert,face,inverse)

try, type = inverse.type;   catch, type = 'GS';     end
try, s    = inverse.smooth; catch, s    = 0.6;      end
try, Np   = inverse.Np;     catch, Np   = 256;      end
try, Nr   = inverse.Nr;     catch, Nr   = 16;       end %% requested number of temporal modes, could be changed depending on svd
try, xyz  = inverse.xyz;    catch, xyz  = [0 0 0];  end
try, rad  = inverse.rad;    catch, rad  = 128;      end
try, hpf  = inverse.hpf;    catch, hpf  = 48;       end %% need to one day put these the correct way round
try, lpf  = inverse.lpf;    catch, lpf  = 0;        end
try, sdv  = inverse.sdv;    catch, sdv  = 4;        end
try, Han  = inverse.Han;    catch, Han  = 1;        end
try, woi  = inverse.woi;    catch, woi  = [];       end
try, Nm   = inverse.Nm;     catch, Nm   = [];       end
try, Nt   = inverse.Nt;     catch, Nt   = [];       end %% fixed number of temporal modes
try, Ip   = inverse.Ip;     catch, Ip   = [];       end
try, QE    = inverse.QE;     catch,  QE=1;          end         %  empty room noise measurement
try, Qe0   = inverse.Qe0;     catch, Qe0   = exp(-5);       end  %% set noise floor at 1/100th signal power i.e. assume amplitude SNR of 10
try, inverse.A;     catch, inverse.A   = [];       end %% orthogonal channel modes

try, SHUFFLELEADS=inverse.SHUFFLELEADS;catch, SHUFFLELEADS=0;end; %% ONLY FOR TESTING - destroyes correspondence between geometry and data

% defaults
%--------------------------------------------------------------------------
type = inverse.type;    % Type of inversion scheme


% get specified modalities to invert (default to all)

%--------------------------------------------------------------------------
Ns = size(L,2);

% A=L';
% UL=A*L;
%Spatial projector
[U,ss,vv]    = spm_svd((L*L'),exp(-16));
A     = U';                 % spatial projector A
UL    = A*L;

AY = A*Y;
AQeA   = A*QE*A';			% Note that here it is A*A'
Qe{1}  = AQeA/(trace(AQeA)); % it means IID noise in virtual sensor space
Q0          = Qe0*trace(Y*Y')*Qe{1}; %% fixed (min) level of sensor space variance


if strcmp(inverse.type_cov,'G')
    AYYAl = AY*AY';
    %Day = squareform(pdist(AY,'mahalanobis'));
    Day = squareform(pdist(AY));
    s0 = median(squareform(Day));
    AYYA = exp(-Day.^2/(2*s0^2)).*AYYAl;
else
    AYYA = AY*AY'; %lineal
end


switch(type)
    case {'MSP','GS','ARD'}
        % create MSP spatial basis set in source space
        %------------------------------------------------------------------
        Qp    = {};
        LQpL  = {};
        Ip    = ceil([1:Np]*Ns/Np);
        for i = 1:Np
            
            % left hemisphere
            %--------------------------------------------------------------
            q               = QG(:,Ip(i));
            Qp{end + 1}.q   = q;
            LQpL{end + 1}.q = UL*q;
            
            % right hemisphere
            %--------------------------------------------------------------
            [d j] = min(sum([vert(:,1) + vert(Ip(i),1), ...
                vert(:,2) - vert(Ip(i),2), ...
                vert(:,3) - vert(Ip(i),3)].^2,2));
            q               = QG(:,j);
            Qp{end + 1}.q   = q;
            LQpL{end + 1}.q = UL*q;
        end
        
    case {'EBB'}
        % create SMOOTH beamforming prior.
        disp('NB smooth EBB algorithm !');
        %------------------------------------------------------------------
        InvCov = spm_inv(AYYA);
        allsource = sparse(Ns,1);
        Sourcepower = sparse(Ns,1);
        for bk = 1:Ns
            q               = QG(:,bk);
            
            smthlead = UL*q;     %% THIS IS WHERE THE SMOOTHNESS GETS ADDED
            normpower = 1/(smthlead'*smthlead);
            Sourcepower(bk) = 1/(smthlead'*InvCov*smthlead);
            allsource(bk) = Sourcepower(bk)./normpower;
        end
        allsource = allsource/max(allsource);   % Normalise
        
        Qp{1} = diag(allsource);
        LQpL{1} = UL*diag(allsource)*UL';
        
    case {'LOR','COH'}
        % create minimum norm prior
        %------------------------------------------------------------------
        Qp{1}   = speye(Ns,Ns);
        LQpL{1} = UL*UL';
        
        % add smoothness component in source space
        %------------------------------------------------------------------
        Qp{2}   = QG;
        LQpL{2} = UL*Qp{2}*UL';
        
    case {'IID','MMN'}
        % create minimum norm prior
        %------------------------------------------------------------------
        Qp{1}   = speye(Ns,Ns);
        LQpL{1} = UL*UL';
end

QP     = {};
LQP    = {};
LQPL   = {};


switch(type)
    
    case {'MSP','GS','EBBgs'}
        % Greedy search over MSPs
        %------------------------------------------------------------------
        Np    = length(Qp);
        Q     = zeros(Ns,Np); %% NB SETTING UP A NEW Q HERE
        for i = 1:Np
            Q(:,i) = Qp{i}.q;
        end
        Q = sparse(Q);
        
        % Multivariate Bayes (Here is performed the inversion)
        %------------------------------------------------------------------
        
        MVB   = spm_mvb(AY,UL,[],Q,Qe,16); %% Qe is identity with unit trace
        
        % Accumulate empirical priors (New set of patches for the second inversion)
        %------------------------------------------------------------------
        % MVB.cp provides the final weights of the hyperparameters
        Qcp           = Q*MVB.cp;
        QP{end + 1}   = sum(Qcp.*Q,2);
        LQP{end + 1}  = (UL*Qcp)*Q';
        LQPL{end + 1} = LQP{end}*UL';
        
end

switch(type)
    
    case {'IID','MMN','LOR','COH','EBB'}
        
        % or ReML - ARD (Here is performed the inversion)
        %------------------------------------------------------------------
        
        
        
        [Cy,h,Ph,F] = spm_reml_sc(AYYA,[],[Qe LQpL],1,-4,16,Q0);
        
        % Spatial priors (QP)
        %------------------------------------------------------------------
        % h provides the final weights of the hyperparameters
        Ne    = length(Qe);
        Np    = length(Qp);
        
        hp    = h(Ne + (1:Np));
        qp    = sparse(0);
        for i = 1:Np
            qp = qp + hp(i)*Qp{i};
        end
        
        % Accumulate empirical priors (New set of patches for the second inversion)
        %------------------------------------------------------------------
        QP{end + 1}   = diag(qp);
        LQP{end + 1}  = UL*qp;
        LQPL{end + 1} = LQP{end}*UL';
end

%==========================================================================
% Step 2: Re-estimate for each subject separately (fusing all modalities)
%==========================================================================

fprintf('Inverting subject 1\n')



% re-do ReML (with informative hyperpriors)
%----------------------------------------------------------------------
Np    = length(LQPL);       % Final number of priors
Ne    = length(Qe);         % Sensor noise prior


Q     = [{Q0} LQPL]; %% sensor corvariance prior:  Qe is identity with unit trace, LQPL is in the units of data

% if rank(AYYA)~=size(A,1),
%     rank(AYYA)
%     size(AYYA,1)
%     warning('AYYA IS RANK DEFICIENT');
% end;



[Cy,h,Ph,F]= spm_reml_sc(AYYA,[],Q,1,-4,16,Q0);



%% recalculate F here

Cp    = sparse(0);
LCp   = sparse(0);
hp    = h(Ne + (1:Np));
for j = 1:Np
    Cp  =  Cp + hp(j)*QP{j};
    LCp = LCp + hp(j)*LQP{j};
end

% MAP estimates of instantaneous sources
%======================================================================
% This is equivalent to M = Cp*UL'*inv(Qe + UL*Cp*UL'))
% with Cp the posterior source covariance (with optimal h values)
M     = LCp'/Cy;

J = M*AY;
%J = J*V';

function [P,x,S,U,W,da] = spm_MDP_select(MDP,varargin)
% aaction selection using active inference
% FORMAT [P,Q,S,U,W] = spm_MDP_select(MDP)
%
% MDP.T           - process depth (the horizon)
% MDP.K           - memory  depth (default 0)
% MDP.N           - number of variational iterations (default 4)
% MDP.S(N,1)      - true initial state
%
% MDP.A(N,N)      - Likelihood of outcomes given hidden states
% MDP.B{M}(N,N)   - transition probabilities among hidden states (priors)
% MDP.C(N,1)      - terminal cost probabilities (prior N over hidden states)
% MDP.D(N,1)      - initial prior probabilities (concentration parameters)
%
% optional:
% MDP.s(1 x T)    -  vector of true states  - for deterministic solutions
% MDP.o(1 x T)    -  vector of observations - for deterministic solutions
% MDP.a(1 x T)    -  vector of action       - for deterministic solutions
% MDP.w(1 x T)    -  vector of precisions   - for deterministic solutions
%
% MDP.G{M}(N,N)   - transition probabilities used to generate outcomes
%                   (default: the prior transition probabilities)
% MDP.B{T,M}(N,N) - transition probabilities for each time point
% MDP.G{T,M}(N,N) - transition probabilities for each time point
%                   (default: MDP.B{T,M} = MDP.B{M})
%
% MDP.plot        - switch to suppress graphics: (default: 0)
%
% produces:
%
% P(M,T)   - probability of emitting action 1,...,M at time 1,...,T
% Q(N,T)   - an array of conditional (posterior) expectations over N hidden
%            states and time 1,...,T
% S(N,T)   - a sparse matrix of ones, encoding the states at time 1,...,T
% U(M,T)   - a sparse matrix of ones, encoding the action at time 1,...,T
% W(1,T)   - posterior expectations of precision
%
% This routine provides solutions of active inference (minimisation of
% variational free energy)using a generative model based upon a Markov
% decision process. This model and inference scheme is formulated
% in discrete space and time. This means that the generative model (and
% process) are  finite state  machines or hidden Markov models whose
% dynamics are given by transition probabilities among states. For
% simplicity, we assume an isomorphism between hidden states and outcomes,
% where the likelihood corresponds to a particular outcome conditioned upon
% hidden states. Similarly, for simplicity, this routine assumes that action
% and hidden controls are isomorphic. If the dynamics of transition
% probabilities of the true process are not provided, this routine will use
% the equivalent probabilities from the generative model.
 
% This particular scheme is designed for situations in which a particular
% action is to be selected at a particular time. The first control state
% is considered the baseline or default behaviour. Subsequent control
% states are assumed to be selected under the constraint that only one
% action can be emitted once. This provides a particular simplification for
% the generative model, allowing a fairly exhaustive model of potential
% outcomes � eschewing a mean field approximation over successive
% control states. In brief, the agent simply represents the current state
% and states in the immediate and distant future.
%
% The transition probabilities are a cell array of probability transition
% matrices corresponding to each (discrete) the level of the control state.
%
% Mote that the conditional expectations are functions of time but also
% contain expectations about fictive states over time at each time point.
% To create time dependent transition probabilities, one can specify a
% function in place of the transition probabilities under different levels
% of control.
%
% Partially observed Markov decision processes can be modelled by
% specifying a likelihood (as part of a generative model) and absorbing any
% probabilistic mapping between (isomorphic) hidden states and outcomes
% into the transition probabilities G.
%
% See also spm_MDP which uses multiple future states and a the mean
% field, approximation for control states � but allows for different actions
% at all times (as in control problems).
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_MDP_select.m 5280 2013-02-24 22:10:36Z karl $
 
% set up and preliminaries
%==========================================================================
 
% options and precision defaults
%--------------------------------------------------------------------------
try, PLOT = MDP.plot; catch, PLOT = 0; end
try, K    = MDP.K;    catch, K = 0;    end
try, N    = MDP.N;    catch, N = 4;    end
 
if PLOT
    if ishandle(PLOT)
        figure(PLOT); clf
        PLOT = 2;
    else
        spm_figure('GetWin','MDP'); clf
    end
end
 
% generative model and initial states
%--------------------------------------------------------------------------
T     = MDP.T;                     % process depth (the horizon)
Ns    = size(MDP.B{1},1);          % number of hidden states
Nb    = size(MDP.B,1);             % number of time-dependent probabilities
Nu    = size(MDP.B,2);             % number of hidden controls
p0    = eps;                       % smallest probability
 
% likelihood model (for a partially observed MDP implicit in G)
%--------------------------------------------------------------------------
try
    A  = MDP.A + p0;
    No = size(MDP.A,1);           % number of outcomes
catch
    A  = speye(Ns,Ns) + p0;
    No = Ns;
end
A     = A*diag(1./sum(A));
lnA   = log(A);
 
 
% transition probabilities (priors)
%--------------------------------------------------------------------------
for i = 1:T
    for j = 1:Nu
        if i == 1 || Nb == T
            B{i,j}   = MDP.B{i,j} + p0;
            B{i,j}   = B{i,j}*diag(1./sum(B{i,j}));
            lnB{i,j} = log(B{i,j});
        else
            B{i,j}   = B{1,j};
            lnB{i,j} = lnB{1,j};
        end
    end
end
 
% terminal cost probabilities (priors)
%--------------------------------------------------------------------------
C     = spm_vec(MDP.C) + p0;
C     = C/sum(C);
lnC   = log(C);

% oncentration parameters of Dirichlet prior over initial state
%--------------------------------------------------------------------------
try
    D  = spm_vec(MDP.D) + 1;
catch
    D  = 2;
end

% generative process (assume the true process is the same as the model)
%--------------------------------------------------------------------------
try
    G = MDP.G;
catch
    G = MDP.B;
end
Ng    = size(G,1);
for i = 1:T
    for j = 1:Nu
        if i == 1 || Ng == T
            G{i,j} = G{i,j} + p0;
            G{i,j} = G{i,j}*diag(1./sum(G{i,j}));
        else
            G{i,j} = G{1,j};
        end
    end
end
 
% initial states and outcomes
%--------------------------------------------------------------------------
[p q]  = max(A*MDP.S(:,1));        % initial outcome (index)
s      = find( MDP.S(:,1));        % initial state (index)
a      = sparse(1,1,1,1,T);        % action (index)
o      = sparse(1,1,q,1,T);        % observations (index)
S      = sparse(s,1,1,Ns,T);       % states sampled (1 in K vector)
O      = sparse(q,1,1,No,T);       % states observed (1 in K vector)
U      = sparse(1,1,1,Nu,T - 1);   % action selected (1 in K vector)
P      = sparse(Nu,T - 1);         % posterior beliefs about control

% hyperpriors
%--------------------------------------------------------------------------
alpha  = 8;
beta   = alpha/8;
W      = alpha/beta;
 
% solve
%==========================================================================
 
% sufficient statistics of hidden states (past, current and last)
%--------------------------------------------------------------------------
x      = zeros(Ns,    T);
u      = zeros(Nu - 1,T);
d      = D/sum(D);
x(:,1) = d;
x(:,T) = C;
da     = [];
for t  = 1:(T - 1)
    
    
    % conditional KL divergence (from C) under available action
    %----------------------------------------------------------------------
    for j = 1:(Nu - 1)
        for k = t:T
            
            % composoiton of future states
            %--------------------------------------------------------------
            p = 1;
            for i = t:T
                if i == k
                    p = B{i,j + 1}*p;
                else
                    p = B{i,1}*p;
                end
            end
            
            % divergence
            %--------------------------------------------------------------
            if nargin > 1
                KL(j,k,:) = -lnC'*p;
            else
                KL(j,k,:) = sum(p.*log(p)) - lnC'*p;
            end
            
        end
    end
    
    
    % Variational iterations (assuming precise inference about past action)
    %----------------------------------------------------------------------
    for i  = 1:N
        
        % past states
        %------------------------------------------------------------------
        for k = max(t - K,1):(t - 1)
            p = max(k - 1,1);
            v = lnA(o(k),:)';
            v = v + lnB{k,a(k)}'*x(:,k + 1);
            v = v + lnB{p,a(p)} *x(:,p    );
            x(:,k) = spm_softmax(v);
        end
        
        % intial state (with a Dirichlet prior)
        %------------------------------------------------------------------
        if t == 1
            v  = x(:,t) + D - 1;
            d  = v/sum(v);
        end       

        % present state
        %------------------------------------------------------------------
        k = t:T;
        p = max(t - 1,1);
        v = lnA(o(t),:)' + lnB{p,a(p)}*x(:,p);
        if t == 1
            v = v + log(d);
        end
        for j = 1:(Nu - 1)
            v = v - W(t)*(u(j,k)*squeeze(KL(j,k,:)))';
        end
        x(:,t) = spm_softmax(v);
        
        
        % precision (W)
        %------------------------------------------------------------------
        try
            W(t)  = MDP.w(t);
        catch
            v     = beta;
            for j = 1:(Nu - 1)
                v = v + u(j,k)*squeeze(KL(j,k,:))*x(:,t);
            end
            W(t)  = alpha/v;
        end
        da(end + 1) = W(t);
        
        
        % policy (u)
        %------------------------------------------------------------------
        v     = 0;
        for j = 1:Ns
            v = v - W(t)*squeeze(KL(:,k,j))*x(j,t);
        end
        u(:,k) = spm_unvec(spm_softmax(spm_vec(v)),v);
        u(:,k) = u(:,k)*(1 - any(a > 1));
 
        
        % graphics to inspect update scheme
        %==================================================================
        
        % save posterior expectations (control)
        %------------------------------------------------------------------
        P  = [1 - sum(u,1); u];
        if PLOT > 2
            spm_plot_states(x,P)
        end
        
    end
    
    
    % graphics to inspect posterior beliefs about hidden states
    %======================================================================
    if PLOT > 1
        spm_plot_states(x,P)
    end
    
    
    % sampling of next state (outcome)
    %======================================================================
    
    % next action (the action that minimises expected free energy)
    %----------------------------------------------------------------------
    try
        a(t) = MDP.a(t);
    catch
        a(t) = find(rand < cumsum(P(:,t)),1);
    end
    
    % next sampled state
    %----------------------------------------------------------------------
    try
        s(t + 1) = MDP.s(t);
    catch
        s(t + 1) = find(rand < cumsum(G{t,a(t)}(:,s(t))),1);
    end
    
    % next obsverved state
    %----------------------------------------------------------------------
    try
        o(t + 1) = MDP.o(t);
    catch
        o(t + 1) = find(rand < cumsum(A(:,s(t + 1))),1);
    end
    
    
    % save action and state sampled
    %----------------------------------------------------------------------
    W(1,t + 1)        = W(t);
    U(a(t)    ,t)     = 1;
    O(o(t + 1),t + 1) = 1;
    S(s(t + 1),t + 1) = 1;

   
    
    % plot
    %======================================================================
    if PLOT > 0
        
        % true state (outcome)
        %------------------------------------------------------------------
        subplot(4,2,5)
        if size(S,1) > 512
            spy(S,16)
        else
            imagesc(1 - S)
        end
        title('True states','FontSize',16)
        ylabel('State','FontSize',12)
        
        % sample (observation)
        %------------------------------------------------------------------
        subplot(4,2,7)
        if size(O,1) > 512
            spy(O,16)
        else
            imagesc(1 - O)
        end
        title('Observed states','FontSize',16)
        xlabel('Time','FontSize',12)
        ylabel('State','FontSize',12)
        
        
        % action sampled (selected)
        %------------------------------------------------------------------
        subplot(4,2,6)
        if size(U,1) > 512
            spy(U,16)
        else
            imagesc(1 - U)
        end
        title('Selected action','FontSize',16)
        ylabel('Action','FontSize',12)
        
        % expected action
        %------------------------------------------------------------------
        subplot(4,2,8)
        plot((1:length(da))/N,da)
        title('Expected precision','FontSize',16)
        xlabel('Time','FontSize',12)
        ylabel('Precision','FontSize',12)
        spm_axis tight
        drawnow
        
        
    end
    
end


% remove last posterior
%--------------------------------------------------------------------------
P = P(:,1:(T - 1));
 
 
function spm_plot_states(x,u)
 
% posterior beliefs about hidden states
%--------------------------------------------------------------------------
subplot(2,2,1)
imagesc(1 - x)
axis square
title('Inferred states','FontSize',16)
xlabel('Time','FontSize',12)
ylabel('Hidden state','FontSize',12)
 
 
% posterior beliefs about control states
%==========================================================================
subplot(2,2,2), hold on
 
% make previous plots dotted lines
%--------------------------------------------------------------------------
h     = get(gca,'Children');
for i = 1:length(h)
    set(h(i),'LineStyle',':');
end
plot(u')
title('Inferred policy','FontSize',16)
xlabel('Time','FontSize',12)
ylabel('Control state','FontSize',12)
spm_axis tight square
legend({'Stay','Shift'}), hold off
drawnow

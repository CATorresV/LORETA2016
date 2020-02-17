
function [yn,x,tint]=eeg_dynamic_simulation(N,inter,SNR,Nd,Nas)
%forward eeg simulation
%yn:noisy output
%x:clean states
%tint:time interval
%Eduardo Giraldo
%egiraldos@utp.edu.co

%carga matrices de campo guia y de vecinos
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
tint = linspace(inter(1),inter(2),N);
%tamaño matriz de campo
Dsize = size(M);

%selecciona fuente aleatoria en la que se quiere simular la actividad
nf=randi([1 Dsize(2)],1,Nas);
%vector de dispersion para simular fuente
v=QG(:,nf);
%actividad en la fuente
ac=0.5*rand(Nas, length(tint));
%parametros
tau=20;
w1=[1.0628;-0.42857;0.008;0.000143;-0.000286];
w2=[1.3;-1;0.008;0.000143;-0.000286];
%inicializacion
x1=zeros(Dsize(2),1);
x2=zeros(Dsize(2),1);
for k=2:N,
    %Nonlinear Model State Matrix
    Gk=[x1 x2 zeros(Dsize(2),1) x1.^2 x1.^3];
    if k>tau,
        Gk(:,3)=x(:,k-tau);
    end
    if k<125
        w=w1;
    else
        w=w2;
    end
    x(:,k)=Gk*w+v*ac(:,k);
    x1=x(:,k);
    if k>=2,
        x2=x(:,k-1);
    end
end
%multiplica por matriz de campo
y=M*x;

yn=y;
for in=1:size(y,1),
    yn(in,:)=awgn(y(in,:),SNR,'measured');
end


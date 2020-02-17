
function [y]=filtro(x,LowPassFreq,HighPassFreq,fs,ord)

%************INPUT DEFINITION**********************************************
%function [y]=filtro(x,LowPassFreq,HighPassFreq,fs);
%x: segnale da filtrare
%LowPassFreq: frequenza passa-basso
%HigPassFreq: frequenza passa-alto
%fs: frequenza di campionamento
%ord: filter order
%**************************************************************************

% x=sin([0.01:0.01:50]')+wgn(5000,1,1);
  
% ord=400;%filter order
% x=SegProva16;
%trasformo le frequenze in hertz in freq normalizzate
hpf=HighPassFreq/(fs/2);
lpf=LowPassFreq/(fs/2);

if hpf==0
    b = fir1(ord,lpf);
elseif lpf==0
    b = fir1(ord,hpf,'high');
else
    b = fir1(ord,[lpf hpf]);
end%if hpf==0

% hd = dfilt.dffir(b);
% freqz(hd);
y=filtfilt(b,1,x);

return
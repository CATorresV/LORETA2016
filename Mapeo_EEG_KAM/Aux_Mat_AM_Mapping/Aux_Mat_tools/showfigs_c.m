%SHOWFIGS Show all figures on the screen
%
%  SHOWFIGS(K)
%
% Use K figures on a row

function showfigs_c(k)

h = sort(get(0,'children'));
n = length(h);
if nargin == 0
	k = ceil(sqrt(n));
end
s = 0.95/k;   % screen stitch
r = 0.93;     % image size reduction
t = 0.055;    % top gap
b = 0.005;    % border gap
fig = 0;
for i=1:k
	for j=1:k
		fig = fig+1;
		if fig > n, break; end
		set(h(fig),'units','normalized','position',[(j-1)*s+b,(1-t)-i*s,s*r,s*r]);
		figure(h(fig));
	end
end
for j=n:-1:1, figure(h(j)); end
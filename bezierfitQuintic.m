% fit one vector (quartic: 4th order)
% t = [x1 x2 x3 x4 ...]
% output: [p2; p3]
function [yy ctrl] = bezierfitQuintic(t,y, tt) % 5th degree

N = length(t);

if (size(y,2)~=N)
    error ('y size should be 1xN (1)')
end
if (size(y,1)~=1)
    error ('y size should be 1xN (2)')
end
if (size(y)~=size(t))
    error ('size should be 1xN (3)')
end

if t(1)~=0
    t=t-t(1);
end
t = t / max(t);

T = [t'.^5  t'.^4  t'.^3  t'.^2   t'  ones(N,1)];
M = [5 -10 10 -5; -20 30 -20 5; 30 -30 10 0; -20 10 0 0 ; 5 0 0 0 ; 0 0 0 0];
Mc= [y(N)-y(1); 5*y(1); -10*y(1); 10*y(1); -5*y(1); y(1)];
TM = T*M;
ctrl = inv(TM'*TM) * TM' * (y'-T*Mc);

yy = TM*ctrl + T*Mc;
end
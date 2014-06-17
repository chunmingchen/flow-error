% fit one vector (quartic: 4th order)
% t = [x1 x2 x3 x4 ...]
% output: [p2; p3]
function [yy ctrl] = bezierfitquartic(t,y, tt)

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

T = [t'.^4  t'.^3  t'.^2   t'  ones(N,1)];
M = [-4 6 -4;  12 -12 4;  -12 6 0;  4 0 0;  0 0 0];
Mc= [y(N)+y(1); -4*y(1); 6*y(1); -4*y(1); y(1)];
TM = T*M;
ctrl = inv(TM'*TM) * TM' * (y'-T*Mc);

yy = TM*ctrl + T*Mc;
end
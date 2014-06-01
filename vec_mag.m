% [x1 x2 x3]
% [y1 y2 y3]
% output : sqrt([x1^2+y1^2...])
function y = vec_mag(x)
    y = sqrt(x(1,:).^2 + x(2,:).^2 + x(3,:).^2);
end
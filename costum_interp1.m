function vq = costum_interp1(x, v, xq, method)
% Extending the interp1 function so that it can handle a group of
% 1-dimensional data.


s = size(v);
vq = zeros([s(1), length(xq)]);

for k = 1 : s(1)
    vq(k, :) = interp1(x, v(k,:), xq, method);
end

end
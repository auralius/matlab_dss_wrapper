function d = zoh_upsample_trim(d, r, n)

d  = conv(upsample(d, r), ones(r,1)); % upsample by r

if nargin > 2
    d = d(1:n); % trim, from element-1 to element-n
end

end
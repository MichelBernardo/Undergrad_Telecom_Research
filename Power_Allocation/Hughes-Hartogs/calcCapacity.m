function [y] = calcCapacity(x)

y = 0;

for iter = 1:length(x)
   
    snr = x(iter)/(sum(x.^2)-x(iter));
    
    y = y + log2(1+snr);
    
end

end
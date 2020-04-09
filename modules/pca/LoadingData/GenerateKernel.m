% Generates a kernel around a given function - TPM 24/02/2020

%k_function = @(distance,r) ((1 - (distance)./r).^2).^2;
%r=3;

function kernel = GenerateKernel(k_function,r)
    
    kernel=zeros(r,r);
    % get distance kernel
    for i=1:size(kernel,1)
        for j=1:size(kernel,2)     
            d(i,j)=(size(kernel,1)/2-i+0.5)^2+(size(kernel,2)/2-j+0.5)^2;
            
            if d(i,j) <=r
                kernel(i,j)=k_function(d(i,j),r);
            end
            
        end
    end
    
    %normalise the kernel (needs to sum to one as patterns are
    %added)
    kernel=kernel./(sum(kernel(:)));

end
function [R,C,reg] = DFT_reg(stack,target,upscale)
    N = size(stack,3);
    reg = zeros(size(stack));
    
    R = zeros(N,1);
    C = zeros(N,1);
    
    %forward
    for i = 1:N
        source = stack(:,:,i);
        S = dftregistrationAlex(fft2(target),fft2(source),upscale);
        R(i) = S(1);
        C(i) = S(2);
        reg(:,:,i) = imtranslate(source,[S(2),S(1)]);
    end
  
end
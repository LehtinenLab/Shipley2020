function [R,C,reg] = DFT_rect(vol, start, upscale)
    Nz = size(vol,3);
    reg = zeros(size(vol));
    target = vol(:,:,start);
    
    R = zeros(Nz,1);
    C = zeros(Nz,1);
    
    %forward
    for i = start:Nz
        source = vol(:,:,i);
        S = dftregistrationAlex(fft2(target),fft2(source),upscale);
        R(i) = S(1);
        C(i) = S(2);
        target = imtranslate(source,[S(2),S(1)]);
        reg(:,:,i) = target;
    end
    
    %backwards
    target = vol(:,:,start);
    for i = flip(1:start)
        source = vol(:,:,i);
        S = dftregistrationAlex(fft2(target),fft2(source),upscale);
        R(i) = S(1);
        C(i) = S(2);
        target = imtranslate(source,[S(2),S(1)]);
        reg(:,:,i) = target;
    end
end
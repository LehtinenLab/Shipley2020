function [output] = dftregistration3D(buf1ft,buf2ft,usfac)

% usfac = 1;

[nr,nc,np]=size(buf2ft);
Nr = ifftshift(-fix(nr/2):ceil(nr/2)-1);
Nc = ifftshift(-fix(nc/2):ceil(nc/2)-1);
Np = ifftshift(-fix(np/2):ceil(np/2)-1);

CC = ifftn(FTpad(buf1ft.*conj(buf2ft),[usfac*nr,usfac*nc,usfac*np]));
CCabs = abs(CC);    
[~,I] = max(CCabs(:));
[row_shift,col_shift,pl_shift] = ind2sub(size(CCabs),I);
% Now change shifts so that they represent relative shifts and not indices

Nr2 = ifftshift(-fix(nr*usfac/2):ceil(nr*usfac/2)-1);
Nc2 = ifftshift(-fix(nc*usfac/2):ceil(nc*usfac/2)-1);
Np2 = ifftshift(-fix(np*usfac/2):ceil(np*usfac/2)-1);
    
row_shift = Nr2(row_shift)/usfac;
col_shift = Nc2(col_shift)/usfac;
pl_shift = Np2(pl_shift)/usfac;
    
    
% If its only one row or column the shift along that dimension has no
% effect. Set to zero.
if nr == 1,
    row_shift = 0;
end
if nc == 1,
    col_shift = 0;
end  

output=[row_shift,col_shift,pl_shift];
return

function [ imFTout ] = FTpad(imFT,outsize)
% imFTout = FTpad(imFT,outsize)
% Pads or crops the Fourier transform to the desired ouput size. Taking 
% care that the zero frequency is put in the correct place for the output
% for subsequent FT or IFT. Can be used for Fourier transform based
% interpolation, i.e. dirichlet kernel interpolation. 
%
%   Inputs
% imFT      - Input complex array with DC in [1,1]
% outsize   - Output size of array [ny nx] 
%
%   Outputs
% imout   - Output complex image with DC in [1,1]
% Manuel Guizar - 2014.06.02

Nout = outsize;
Nin = size(imFT);
imFT = fftshift(imFT);
center = floor(size(imFT)/2)+1;

imFTout = zeros(outsize);
centerout = floor(size(imFTout)/2)+1;

% imout(centerout(1)+[1:Nin(1)]-center(1),centerout(2)+[1:Nin(2)]-center(2)) ...
%     = imFT;
cenout_cen = centerout - center;

imFTout(max(cenout_cen(1)+1,1):min(cenout_cen(1)+Nin(1),Nout(1)),max(cenout_cen(2)+1,1):min(cenout_cen(2)+Nin(2),Nout(2)),...
    max(cenout_cen(3)+1,1):min(cenout_cen(3)+Nin(3),Nout(3))) ...
    = imFT(max(-cenout_cen(1)+1,1):min(-cenout_cen(1)+Nout(1),Nin(1)),max(-cenout_cen(2)+1,1):min(-cenout_cen(2)+Nout(2),Nin(2)),...
    max(-cenout_cen(3)+1,1):min(-cenout_cen(3)+Nout(3),Nin(3)));

imFTout = ifftshift(imFTout)*Nout(1)*Nout(2)*Nout(3)/(Nin(1)*Nin(2)*Nin(3));
return
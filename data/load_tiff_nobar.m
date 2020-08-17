function img = load_tiff(fpath)
%LOAD_TIFF loads a tiff image from FPATH and outputs a double precison
%array of size(height,width,frames)
%   MJLM 2012-10-08: created

info = imfinfo(fpath,'tiff');

frames = length(info);
width  = info(1).Width;
height = info(1).Height;

img = zeros(height,width,frames);


for i=1:frames
	img(:,:,i) = imread(fpath,'tiff',i,'Info',info);
end
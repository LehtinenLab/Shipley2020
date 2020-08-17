function data = aligned(path, regpath, k, N, pmt, optolevel)
%ALIGNED Aligns a chunk of a file based on regpath
%   Input:
%       path - path to .sbx file, .sbx automatically appended if not there
%       startframe - first frame to read, 1-indexed
%       nframes - number of frames to read
%       tform - a cell array of affine2d transforms to apply first
%       dft - a vector of size(nframes, 4) of dft registration to apply
%       pmt - which color to read, 0 if only one color, 1 if two colors
%           and red
%       [removeedges] - remove the edges before returning, T/F
%   Output:
%       data - an array of size(height, width, min(nframes, max possible))
%           of data that has been registered with affine and dft transforms

    if nargin < 5, pmt = 1; end
    if nargin < 6, optolevel = []; end
    
    % Load image data and registration data
    data = pipe.imread(path, k, N, pmt, optolevel);
    info = pipe.metadata(path);
    reg = load(regpath, '-mat');
    
    if length(size(data)) == 3
        data = reshape(data, 1, size(data, 1), size(data, 2), size(data, 3));
    end
    
    % Get the positions to return, accounting for optotune levels
    pos = k:info.nframes;
    if ~isempty(optolevel) && info.optotune_used
%         pos = optotune_level:length(info.otwave):info.max_idx;
        pos = optolevel:length(info.otwave):info.nframes;
        pos = pos(pos >= k);
%         pos = pos(pos <= k + size(data, 4));
        pos = pos(1:size(data,4));
    end
    
    if isfield(reg, 'tform')
        tform = reg.tform(pos);
        
        for c = 1:size(data, 1)
            for j = 1:size(data, 4)
                disp(j);
                slice = squeeze(data(c, :, :, j));
                data(c, :, :, j) = imwarp(slice, tform{j}, 'OutputView', ...
                    imref2d(size(slice)));
            end
        end
    end
    
    if isfield(reg, 'trans') && (~isfield(reg, 'binframes') || reg.binframes > 1)
        cl = class(data);
        
        fr1 = zeros(info.sz(1), info.sz(2));
        [nr, nc] = size(fft2(double(fr1)));
        Nr = ifftshift(-fix(nr/2):ceil(nr/2)-1);
        Nc = ifftshift(-fix(nc/2):ceil(nc/2)-1);
        [Nc, Nr] = meshgrid(Nc, Nr);
        
        trans = reg.trans(pos, :);
        for j = 1:size(data, 4)
            row_shift = trans(j, 3);
            col_shift = trans(j, 4);
            diffphase = trans(j, 2);

            for c = 1:size(data, 1)
                fftslice = fft2(double(squeeze(data(c, :, :, j))));
                frame = fftslice.*exp(1i*2*pi*(-row_shift*Nr/nr - col_shift*Nc/nc));
                frame = frame*exp(1i*diffphase);
                data(c, :, :, j) = cast(abs(ifft2(frame)), cl);
            end
        end
    end
    
%     if size(data, 1) == 1
%         data = reshape(data, size(data, 2), size(data, 3), size(data, 4));
%     end
end


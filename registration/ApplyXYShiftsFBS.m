function[correctedVolume] = ApplyXYShiftsFBS(unshiftedVolume,...
    RowShifts, ColumnShifts)
    
%   APPLYXYSHITS: apply XY shifts per plane with a parfor loop
%   Credits: Alex Fratzl
%
%   Inputs:
%     correctedVolume -- 4D matrix of uint16 or other, dim (x,y,z,t)
%     RowShifts -- 2D matrix of doubles, dim (z,t)
%     ColumnShifts -- 2D matrix of doubles, dim (z,t)
%   Outputs:
%     correctedVolume -- 4D matrix of uint16 or other, dim (x,y,z,t)

    nbplanes = size(unshiftedVolume, 3); % get number of frames before the
    % parfor loop
    correctedVolume = zeros(size(unshiftedVolume));
    for t = 1:size(unshiftedVolume, 4)
        for i = 1:nbplanes
            R = RowShifts(i,t);
            C = ColumnShifts(i,t);
            slice = unshiftedVolume(:,:,i,t);
            correctedVolume(:,:,i,t) = imtranslate(slice, [C, R]); % careful order!
        end
    end

end
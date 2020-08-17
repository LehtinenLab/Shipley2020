function tforms_optotune = CalculateOptotuneWarp(path,refchannel,scale,varargin)
    p = inputParser;
    addOptional(p,'edges',[0,0,0,0]);
    addOptional(p,'regtype','affine'); %options are 'affine', 'rigid', or 'none'
    addOptional(p,'refsize',30);
    addOptional(p,'save',false);
    
    
    parse(p,varargin{:});
    p = p.Results;
    
% CALCULATE OPTOTUNE WARPING
    info = pipe.io.sbxInfo(path);
    fdir = fileparts(path);
    
    Nx = info.sz(1) - p.edges(3) - p.edges(4);
    Ny = info.sz(2) - p.edges(1) - p.edges(2);
    Nz = size(info.otwave,2);
    
    if  strcmp(p.regtype, 'none')
        tforms_optotune = repmat(affine2d(eye(3)),[1,Nz]);
        return;
    end
    
    %load n volumes, and crops the edges
    raw_ref = pipe.imread(path,1,p.refsize*Nz,refchannel,[]);
    raw_ref = raw_ref(p.edges(3)+1:end-p.edges(4),p.edges(1)+1:end-p.edges(2),:);
    raw_ref = reshape(raw_ref,Nx,Ny,Nz,[]);
    raw_ref = imresize(raw_ref,1/scale);

    %make one mean volume from the 30 volume chunk
    mean_raw_ref = squeeze(mean(raw_ref,4));
    % mean_raw_ref = squeeze(median(raw_ref,4));

    if strcmp(p.regtype,'affine')
    %calculate optotune warping transformation from the mean volume
        tforms_optotune = MultiStackReg_Fiji_affine_2(mean_raw_ref,fdir,Nz);
    elseif strcmp(p.regtype, 'rigid')
        tforms_optotune = MultiStackReg_Fiji_rigid(mean_raw_ref,fdir,Nz);
    else
        fprintf('Invalid registration type for optotune correction');
        return;
    end
    
    %save the transformations for later
    if p.save
        save(strcat(fdir,filesep,'tforms_optotune.mat'),'tforms_optotune');
    end
    
end
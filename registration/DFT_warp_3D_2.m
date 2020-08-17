function DFT_warp_3D_2(path, shiftpath, refchannel, scale, Nchunks, tforms_optotune, varargin)
    
    p = inputParser;
    addOptional(p,'edges',[0,0,0,0]);
    addOptional(p,'Nt',[]); %how many volumes to register (if empty, does all)
    addOptional(p,'optotune','true'); %whether to apply the optotune transformation
    addOptional(p,'reftype','median'); %options are 'median' or 'mean' for projecting the reference volume
    addOptional(p,'blurfactor',1); %width of gaussian to blur for DFT reg
    addOptional(p,'keepingfactor',0.95); %what proportion of frame to accountf or with shifts
    addOptional(p,'planescorr',3); %how many planes up and down to search for Z DFT reg
    addOptional(p,'save',true); %save the DFT and optotune shifts
    
    parse(p,varargin{:});
    p = p.Results;

    info = pipe.io.sbxInfo(path);
    fdir = fileparts(path);
    
    Nx = info.sz(1) - p.edges(3) - p.edges(4);
    Ny = info.sz(2) - p.edges(1) - p.edges(2);
    Nz = size(info.otwave,2);
    
    if isempty(p.Nt)
        p.Nt = info.nframes / numel(info.otwave);
    end
        
    chunkframes = Nz*floor(p.Nt / Nchunks);

    RS = cell(1,Nchunks);
    CS = cell(1,Nchunks);
    ZS = cell(1,Nchunks);
    
    H = parfor_progressbar(Nchunks,'DFT registration');

    
    for chunk = 1:Nchunks
        %1) load reference chunk
        raw_chunk = pipe.imread(path,chunkframes*(chunk-1)+1, chunkframes,refchannel,[]);
        raw_chunk = raw_chunk(p.edges(3)+1:end-p.edges(4),p.edges(1)+1:end-p.edges(2),:);
        raw_chunk = reshape(raw_chunk,Nx,Ny,Nz,[]);
        raw_chunk = imresize(raw_chunk,1/scale);
    
        
        %2) apply optotune registration, if desired
        if p.optotune
            unwarped_chunk = ApplyOptotuneWarp(raw_chunk,tforms_optotune);
        else
            unwarp_chunk = raw_chunk;
        end
        
        % rectify each volume with dft
        RS0 = zeros(Nz,size(unwarped_chunk,4));
        CS0 = zeros(Nz,size(unwarped_chunk,4));
        chunk_reg0 = zeros(size(unwarped_chunk));
        for i = 1:size(unwarped_chunk,4)
            [RS0(:,i),CS0(:,i),chunk_reg0(:,:,:,i)] = DFT_rect(unwarped_chunk(:,:,:,i),round(Nz/2),4);
        end
        
        
        
        %3) first round of XY DFT registration
        %make single target volume for entire chunk
        ref1 = defineReference(chunk_reg0,size(chunk_reg0,4),p.reftype);
        %calculate DFT shift between each volume and target vol
        [RS1, CS1] = DetermineXYShiftsFBS(chunk_reg0,p.blurfactor,p.keepingfactor,ref1);
        %apply the shift (need it for subsequent steps)
        chunk_reg1 = ApplyXYShiftsFBS(chunk_reg0,RS1,CS1);
        
        ref2 = defineReference(chunk_reg1,size(chunk_reg0,4),p.reftype);
        
        shifts = zeros(size(chunk_reg1,4),3);
        for j = 1:size(chunk_reg1,4)
            vol = chunk_reg1(:,:,:,j);
            shifts(j,:) = dftregistration3D(fftn(ref2),fftn(vol),2);
        end
        
        RS2 = shifts(:,1)*scale;
        CS2 = shifts(:,2)*scale;
        ZS1 = shifts(:,3);
        
        
        %6) combine Row and Column Shifts from (3), (4), and (5)
        RS(:,chunk) = {RS0*scale+RS1*scale+repmat(RS2',[Nz,1])};
        CS(:,chunk) = {CS0*scale+CS1*scale+repmat(CS2',[Nz,1])};
        ZS(:,chunk) = {ZS1'};
        

        %7) save reference files for stitching later
        ref_all{chunk} = ref2;
        RS0_all{chunk} = RS0;
        CS0_all{chunk} = CS0;
        RS1_all{chunk} = RS1;
        CS1_all{chunk} = CS1;
        RS2_all{chunk} = RS2;
        CS2_all{chunk} = CS2;
        
        %iterate 
        H.iterate(1);
    end
    
    intermediate_shifts.RS0_all = [RS0_all{:}];
    intermediate_shifts.RS1_all = [RS1_all{:}];
    intermediate_shifts.RS2_all = [RS2_all{:}];
    intermediate_shifts.CS0_all = [CS0_all{:}];
    intermediate_shifts.CS1_all = [CS1_all{:}];
    intermediate_shifts.CS2_all = [CS2_all{:}];
    
    
    %Fix intra-chunk discontinuities
    ref_final = ref_all{1};
    
    interchunk_shifts = zeros(Nchunks,3);
    for j = 1:Nchunks
        interchunk_shifts(j,:) = dftregistration3D(fftn(ref_final),fftn(ref_all{j}),2);
    end
    RS_chunk = interchunk_shifts(:,1)*scale;
    CS_chunk = interchunk_shifts(:,2)*scale;
    ZS_chunk = interchunk_shifts(:,3);
    
    %convert local shift correction cells to matrix form
    RS = [RS{:}];
    CS = [CS{:}];
    ZS = [ZS{:}];
    

    %stretch the intra-chunk corrections to apply to every frame
    RS_chunk = imresize(RS_chunk',size(RS),'nearest');
    CS_chunk = imresize(CS_chunk',size(CS),'nearest');
    ZS_chunk = imresize(ZS_chunk',size(ZS),'nearest');
   
    
    %scale the optotune transforms
    tforms_optotune_full = tforms_optotune;
    for i = 1:size(tforms_optotune,2)
        tforms_optotune_full(i).T(3,1:2) = tforms_optotune_full(i).T(3,1:2)*scale;
    end
    
    %save the DFT and optotune transformations
    if p.save
        save(shiftpath,'RS','CS','ZS','RS_chunk','CS_chunk','ZS_chunk','intermediate_shifts','scale','tforms_optotune_full','ref_all','-mat');
    end
    close(H);
end
























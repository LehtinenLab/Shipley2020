function [zproj_reg,R,C,TF] = zproj_reg(k, Nt, pmt, otrange, varargin)

    p = inputParser;
    addOptional(p, 'pathz',[]);
    addOptional(p, 'proj_type', 'mean');
    addOptional(p, 'mtype', '.sbx'); %movie type suffix (other options are .sbxz or .sbxreg)
    addOptional(p, 'reg_suffix', []); %add a registration path
    addOptional(p, 'align', false); %use pipe.align registration
    addOptional(p, 'scale', 4);
    addOptional(p, 'regtype', 'DFT');
    addOptional(p, 'refchan', 1);
    addOptional(p, 'write_unreg', false);
    addOptional(p, 'zproj_raw', []);
    parse(p, varargin{:});
    p = p.Results;

    if isempty(p.zproj_raw)
        zproj = pipe.zproj(p.pathz, k, Nt, pmt, otrange, 'mtype',p.mtype, 'registration',p.align, 'reg_suffix',p.reg_suffix,'proj_type',p.proj_type);
    else
        zproj = p.zproj_raw;
    end
    
    %% do DFT align
    zproj_ref = squeeze(zproj(p.refchan,:,:,:));
    
    raw_ref = imresize(zproj_ref,1/p.scale);
    
    for i = 1:Nt
        slice = raw_ref(:,:,i);
        L = median(slice(:));
        U = prctile(slice(:),99);
        slice = rescale(slice,'InputMin',L,'inputMax',U);
        raw_adj(:,:,i) = slice;
    end
    
    
    target1 = mean(raw_adj(:,:,1:min(Nt,50)),3);
    [R1,C1,reg1] = DFT_reg(raw_adj,target1,p.scale);
    [R2,C2,reg2] = DFT_rect(reg1,round(Nt/2),p.scale);
    target3 = median(reg2,3);
    [R3,C3,reg3] = DFT_reg(reg2,target3,p.scale);
    
    if strcmp(p.regtype,'Affine')
        TF = MultiStackReg_Fiji_affine(reg3,strcat(fileparts(pathz),filesep),size(reg3,3));
        for i = 1:Nt
            TF(i).T(3,1:2) = p.scale * TF(i).T(3,1:2);
        end
    elseif strcmp(p.regtype,'DFT')
        TF = repmat(affine2d(eye(3)),[1,Nt]);
    end
    
    R = (R1+R2+R3)*p.scale;
    C = (C1+C2+C3)*p.scale;
    zproj_reg = zeros(size(zproj));
    w = waitbar(0,'applying transformations...');
    for i = 1:Nt
        for c = 1:size(zproj,1)
            slice = squeeze(zproj(c,:,:,i));
            slice = imtranslate(slice,[C(i),R(i)]);
            slice = imwarp(slice,TF(i),'OutputView',imref2d(size(slice)));
            zproj_reg(c,:,:,i) = slice;
        end
        waitbar(i/Nt);
    end
    close(w);
end
















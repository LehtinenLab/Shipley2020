function write2chanTiff(mov,path)
    if ndims(mov)==4
        Nc = size(mov,1);
        Nx = size(mov,2);
        Ny = size(mov,3);
        Nt = size(mov,4);

        mov = permute(mov,[2,3,1,4]);
        mov = reshape(mov,Nx,Ny,[]);
        Miji(false);
        MIJ.createImage(mov);
        MIJ.run('Stack to Hyperstack...', sprintf('order=xyczt(default) channels=%d slices=%d frames=%d display=Composite',Nc,1,Nt));
        MIJ.run('Save', strcat('Tiff..., path=[',path,']'));
        MIJ.closeAllWindows
        MIJ.exit;
    elseif ndims(mov) == 5
        Nc = size(mov,1);
        Nx = size(mov,2);
        Ny = size(mov,3);
        Nz = size(mov,4);
        Nt = size(mov,5);

        mov = permute(mov,[2,3,1,4,5]);
        mov = reshape(mov,Nx,Ny,[]);
        Miji(false);
        MIJ.createImage(mov);
        MIJ.run('Stack to Hyperstack...', sprintf('order=xyczt(default) channels=%d slices=%d frames=%d display=Composite',Nc,Nz,Nt));
        MIJ.run('Save', strcat('Tiff..., path=[',path,']'));
        MIJ.closeAllWindows
        MIJ.exit;
    else
        disp('Movie dim must be 4 (CXYT) or 5 (CXYZT) for writing');
    end
end
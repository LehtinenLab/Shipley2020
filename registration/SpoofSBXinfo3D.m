function info = SpoofSBXinfo3D(yDim, xDim, zDim, tDim, nchan)

info = struct;
info.resfreq = 7930;
info.postTriggerSamples = 5000;
info.nchannels = 1;
info.abort_bit = 0;
info.scanbox_version = 2;
info.volscan =  1;
info.opto2pow = [];
info.power_depth_link = 0;
info.area_line = true;

%altered
info.sz = [yDim,xDim];
info.height = info.sz(2);
info.width = info.sz(1);
info.nframes = zDim*tDim;
info.nchan = nchan;
info.max_idx = info.nframes-1;
info.recordsPerBuffer = info.width;
info.bytesPerBuffer = info.postTriggerSamples * info.width * 2 * info.nchan;
info.nsamples = info.width * info.height * 2 * info.nchan;
if nchan == 2
    info.channels = 1;
else
    info.channels = 2;
end
info.scanmode = 1;


info.optotune_used = 1;
info.otlevels = zDim;
info.otwave = 1:info.otlevels;
end
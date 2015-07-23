function write_ATL06_h5(D3, h5_file);


if exist(h5_file,'file'); delete(h5_file); end

ff=fieldnames(D3);
for kf=1:length(ff);
    temp=D3.(ff{kf});
    this_field=['/', ff{kf}];
    h5create(h5_file, this_field, size(temp),'ChunkSize', [min(size(temp,1), 1024), 1], 'Datatype','double','Deflate', 9);
    h5write(h5_file, this_field,  temp);
end

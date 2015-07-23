function D3=read_ATL06_h5(h5_file)


I=h5info(h5_file,'/');
for kD=1:length(I.Datasets);
    D3.(I.Datasets(kD).Name)=h5read(h5_file,['/', I.Datasets(kD).Name]);
end

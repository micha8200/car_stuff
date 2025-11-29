function str = file2string(filename)
fid = fopen(filename,'r');
if fid>0
    str = fread(fid, '*char')';
end
fclose(fid);
end
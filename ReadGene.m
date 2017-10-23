function data = ReadGene(fileName)
    fid=fopen(fileName,'rb');
    if(fid>0)
        data = fread(fid,inf,'uint8');
    end
    fclose(fid);
end
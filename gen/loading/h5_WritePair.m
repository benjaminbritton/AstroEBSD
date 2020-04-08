function h5_WritePair(h5_name,h5_gname,dtype,h5_dname,data1)
%H5_WRITEPAIR Write EBSD Data to file

try
    data_name=['/' h5_gname dtype h5_dname];
    
    try 
        h5info(h5_name,data_name);
        disp(['Overwriting ' data_name ]);
    catch
        h5create(h5_name,data_name,size(data1));
    end
    h5write(h5_name,data_name,data1);
catch err
    disp('The data did not write')
    rethrow(err)
end

end


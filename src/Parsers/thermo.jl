function parsethermo(csvfilename::AbstractString)::QPCRDataset
    tempbuf=IOBuffer()
    filelines=split(read(csvfilename,String),"\n")
    #find the end of the header
    linenumber=1
    for line in filelines
        if split(line,",")[1] == "Well"
            #this is the line containing our column names
            break
        end
        linenumber+=1
    end
    filetextnoheader=join(filelines[linenumber:end],"\n")
    write(tempbuf,filetextnoheader)
    seek(tempbuf,0)
    rawdata=CSV.read(tempbuf,DataFrame)
    close(tempbuf)
    data=select(rawdata,"Sample Name"=>:sample,"Target Name"=>:target,"CT"=>ByRow(function(ct)
                                                                                      if ct=="Undetermined"
                                                                                          return Float64(40)
                                                                                      else
                                                                                          return parse(Float64,ct)
                                                                                      end
                                                                                  end)=>:ct)
    return QPCRDataset(data)
end

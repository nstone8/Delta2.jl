#function for testing if a character is whitespace or a quote
function shouldremove(c::Char)::Bool
    quotechars=['"','\'']
    if c in quotechars
        return true
    elseif isspace(c)
        return true
    end
    return false
end

#wizard for performing DeltaCT
function DeltaCT()

end

function DeltaCT(filename::AbstractString)
    
end


DeltaCT(data::QPCRDataset,hkgenes::Vector{String};ntc="NTC",ampthreshold=40)::DeltaCT

DDCT(dct::DeltaCT,refsample)

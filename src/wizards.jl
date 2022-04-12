#function for testing if a character is whitespace or a quote
function escapechars(c::Char)::Bool
    quotechars=['"','\'']
    if c in quotechars
        return true
    elseif isspace(c)
        return true
    end
    return false
end

#wizard for performing DeltaCT
function DeltaCT(qpcr::QPCRDataset)

    #DeltaCT(data::QPCRDataset,hkgenes::Vector{String};ntc="NTC",ampthreshold=40)::DeltaCT
    
end

#wizard for performing DDCT
function DDCT(qpcr::QPCRDataset)

    #DDCT(dct::DeltaCT,refsample)
end

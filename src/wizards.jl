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

function saveprompt(result)::String
    dirprompt="""
              Enter the path to a directory where you would like your results to be saved.
              Tip: Try dragging a folder onto this window
              """
    println(dirprompt)
    savedir=strip(escapechars,readline())
    nameprompt="Enter a filename or press enter to name the file qpcr_results.csv"
    println(nameprompt)
    savefile=nothing
    while true
        savefile=strip(escapechars,readline())
        if (length(savefile)<4) || (savefile[end-3:end]!=".csv")
            savefile*=".csv"
        end
        if isfile(joinpath(savedir,savefile))
            println()
            println("""$(joinpath(savedir,savefile)) already exists.
                        Please choose a different name""")
            
        else
            break             
        end
        
    end
    result |> writeresult(joinpath(savedir,savefile))
end

#wizard for performing DeltaCT
function DeltaCT(qpcr::QPCRDataset,promptsave=true)
    println("select housekeeping genes")
    alltargets=targets(qpcr) |> collect
    println()
    hkgene_indices=MultiSelectMenu(alltargets) |> request |> collect
    println()
    hkgenes=alltargets[hkgene_indices]
    println("was there a no template control on this plate?")
    hasntc=request(RadioMenu(["yes","no"])) == 1 ? true : false
    println()
    ntc=missing
    ampthreshold=40
    if hasntc
        allsamples=samples(qpcr) |> collect
        println("select no template control")
        ntc_index=RadioMenu(allsamples) |> request
        println()
        ntc=allsamples[ntc_index]
        println("set the amplification threshold")
        ampthreshold=41-(string.(40:-1:1) |> RadioMenu |> request)
        println()
    end
  
    println("performing ΔCT with the following settings:")
    println("normalizing expression against targets: "*join(hkgenes,", ", " and "))
    if ismissing(ntc)
        println("not checking for amplification of a no template control")
    else
        println("checking that no targets have a CT value of greater than $ampthreshold for sample $ntc")
    end
    
    dct=DeltaCT(qpcr,hkgenes,ntc=ntc,ampthreshold=ampthreshold)
    println()
    if promptsave
        saveprompt(dct)
    end
    return dct
end

#wizard for performing DDCT
function DDCT(qpcr::QPCRDataset)
    dct=DeltaCT(qpcr,false)
    allsamples=samples(dct) |> collect
    println("select reference sample")
    refindex=RadioMenu(allsamples) |> request
    println()
    refsample=allsamples[refindex]
    ddct=DDCT(dct::DeltaCT,refsample)
    ddct |> saveprompt
    return ddct
end

function askforpath()
    prompt="""
    please enter the path to your qpcr data
    this must be a .csv file containing columns that identify samples, targets and CT values
    tip: try dragging the file onto this window
    """
    println(prompt)
    response=strip(Delta2.escapechars,readline())
    println()
    return response
end

function qpcrdatasetwizard()
    println()
    rawframe=CSV.read(askforpath(),DataFrame)
    allcols=names(rawframe)
    
    println("select column in dataset which contains sample names")
    samp_col_index=RadioMenu(allcols) |> request
    println()

    println("select column in dataset which contains target names")
    target_col_index=RadioMenu(allcols) |> request
    println()

    println("select column in dataset which contains ct values")
    ct_col_index=RadioMenu(allcols) |> request
    println()
    return QPCRDataset(rawframe,allcols[samp_col_index],allcols[target_col_index],allcols[ct_col_index],dropmissing=true)
end
    
#wizard for performing DeltaCT
"""
```julia
DeltaCT()
```
Start a wizard for performing ΔCT on qPCR data contained in a generic .csv file
"""
function DeltaCT()
    qpcrdatasetwizard() |> DeltaCT
end

#wizard for performing DDCT
"""
```julia
DDCT()
```
Start a wizard for performing ΔΔCT on data contained in a generic .csv file
"""
function Delta2.DDCT()
    qpcrdatasetwizard() |> DDCT
end

struct QPCRDataset
    data::AbstractDataFrame
    #add an inner constructor to check that the column labels are right
    function QPCRDataset(data::AbstractDataFrame)::QPCRDataset
        @assert names(data)==["sample","target","ct"] "Incorrect column names, users should not call this constructor directly. Please use a parser function in the Delta2.Parsers submodule"
        @assert eltype(data.ct) <: Real "CT values not numeric. Users should not call this constructor directly. Please use a parser function in the Delta2.Parsers submodule"
        new(data)
    end
end

struct DeltaCT
    data::AbstractDataFrame
    hkgenes::Vector{String}
    function DeltaCT(data::QPCRDataset,hkgenes::Vector{String};ntc="NTC",ampthreshold=40)::DeltaCT
        dataframe=data.data #unwrap our data from the QPCRDataset
        targets=Set(dataframe.target)
        #make sure all the genes listed in hkgenes are represented in targets
        @assert map(hkgenes) do hk
            hk in targets
        end |> all "Not all provided housekeeping genes are present in data"

        #check to see if there was amplification in any of the no template controls
        if !ismissing(ntc)
            ntc_data=dataframe[dataframe.sample.==ntc,:]
            @assert size(ntc_data)[1] > 0 "No data found for no template controls. If there are not any no template controls please pass `ntc=missing`"
            amplified_ntc=ntc_data[ntc_data.ct .< ampthreshold,:]
            if size(amplified_ntc)[1] > 0 #there was amplification in a no template control
                error("Amplification detected in no template controls for the following targets: " * join(amplified_ntc.target," ,"))
            end
        end

        #do deltact
        if !ismissing(ntc)
            ampdata=dataframe[dataframe.sample .!= ntc,:]
        else
            ampdata=copy(dataframe)
        end
        #Take the mean of all our technical replicates
        ampdata=combine(groupby(ampdata,[:sample,:target]),:ct => mean =>:ct)

        #pull out housekeeping data
        hkdata=ampdata[collect(t in hkgenes for t in ampdata.target),:]

        #take the geomean of all our housekeeping genes to get a normalization factor
        dct_norm=combine(
            groupby(hkdata[:,[:sample,:ct]],:sample ),
            :ct => geomean =>:dct_normfactor
        )
        #add our normalization data as a new column to our ampdata
        dct_frame=innerjoin(ampdata,dct_norm,on=:sample)
        #perform the normalization
        dct_frame.dct=dct_frame.ct .- dct_frame.dct_normfactor
        #add a fold change column
        dct_frame.dct_fold=2 .^ (-1 .* dct_frame.dct)
        return new(dct_frame,hkgenes)
    end
end

struct DDCT
    data::AbstractDataFrame
    dct::DeltaCT
    refsample
    function DDCT(dct::DeltaCT,refsample)
        dct_frame=dct.data
        #get the dct values for our reference sample for each target
        ddct_norm=select(dct_frame[dct_frame.sample .== refsample,:], :target, :dct => :ddct_normfactor)

        #add this ddct normalization factor to our existing dct data as a new column
        ddct_frame=innerjoin(dct_frame,ddct_norm,on=:target)

        #perform the normalization
        ddct_frame.ddct=ddct_frame.dct .- ddct_frame.ddct_normfactor

        #add a ddct_fold column
        ddct_frame.ddct_fold = 2 .^ (-1 .* ddct_frame.ddct)
        return new(ddct_frame,dct,refsample)
    end
end

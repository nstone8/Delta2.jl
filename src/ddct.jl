#create an abstract supertype for our analysis results
abstract type DeltaResult end

#function for selecting specific targets and samples frome a DeltaCT.data or DDCT.data DataFrame
function selectdata(data,selectedtargets,selectedsamples)::DataFrame
    targetrows=[r in selectedtargets for r in data.target]
    samplerows=[r in selectedsamples for r in data.sample]
    return data[targetrows .&& samplerows,:] |> copy #wrap in a copy() so I can sleep at night
end

"""
Type representing the results of a qPCR experiment. These objects should be created by a parser
function such as `readpcr`
"""
struct QPCRDataset <: DeltaResult
    data::AbstractDataFrame
    #add an inner constructor to check that the column labels are right
    function QPCRDataset(data::AbstractDataFrame)::QPCRDataset
        @assert names(data)==["sample","target","ct"] "Incorrect column names, users should not call this constructor directly."
        @assert ((eltype(data.ct) <: Real) && (eltype(data.sample) <: AbstractString) && (eltype(data.target) <: AbstractString)) "Column datatypes not correct. Users should not call this constructor directly."
        new(data)
    end
    
    #another inner constructor for indexing
    function QPCRDataset(pcr::QPCRDataset,selectedtargets::Vector{String},selectedsamples::Vector{String})
        st=Set(selectedtargets)
        ss=Set(selectedsamples)
        newdata=selectdata(pcr.data,st,ss)
        new(newdata)
    end
    
end

"""
```julia
QPCRDataset(data, sample_col, target_col, ct_col; [noamp_flag], [dropmissing])
```
Create a QPCRDataset from a .csv file. If the data contains a flag for samples that didn't amplify 
(for instance, an entry of "Undetermined" in the ct column) this flag can be passed to the
`noamp_flag` keyword argument. If the argument dropmissing=true is provided, all
rows with empty values will be dropped after selecting columns.
"""
function QPCRDataset(data::DataFrame, sample_col, target_col, ct_col; noamp_flag=nothing, dropmissing=false)::QPCRDataset
    if dropmissing
        select!(data,sample_col,target_col,ct_col)
        DataFrames.dropmissing!(data)
    end
    selected_data=select(data,sample_col => function(s)
                             #make sure sample names are always Strings
                             String.(s)
                         end => :sample,target_col=>function(t)
                             #make sure targets are also strings
                             String.(t)
                         end => :target,ct_col=>ByRow(function(ct)
                                                          if !isnothing(noamp_flag)
                                                              if ct==noamp_flag
                                                                  return Float64(40)
                                                              end
                                                          end
                                                          if !isa(ct,Real)
                                                              return parse(Float64,ct)
                                                          end
                                                          return ct
                                                      end)=>:ct)
    return QPCRDataset(selected_data)
end

#create an empty readpcr function that other packages implementing parsers can add methods to
function readpcr end


"""
```julia
targets(result)
```
Get the targets present in a `DeltaResult`
"""
targets(dr::DeltaResult) = Set(dr.data.target)

"""
```julia
samples(result)
```
Get the samples present in a `DeltaResult`
"""
samples(dr::DeltaResult) = Set(dr.data.sample)

#make DeltaResults print well in the REPL
function Base.show(io::IO,x::DeltaResult)
    show(io,x.data[:,[1,2,3,end]])
end

#need a more specific option for QPCRDataset
function Base.show(io::IO,x::QPCRDataset)
    show(io,x.data)
end

#support indexing
function Base.getindex(dr::T,targetselector,sampleselector)::T where {T<:DeltaResult}
    st=doselection(targets(dr),targetselector)
    ss=doselection(samples(dr),sampleselector)
    T(dr,st,ss)
end

#allow people to only specify targets
function Base.getindex(dr::T,targetselector)::T where {T<:DeltaResult}
    st=doselection(targets(dr),targetselector)
    ss=samples(dr) |> collect
    T(dr,st,ss)
end

#do selection for Set{String}
function doselection(options::Set{String},selector::Set{String})::Vector{String}
    intersect(options,selector) |> collect
end

#do selection for Vector{String}
function doselection(options::Set{String},selector::Vector{String})::Vector{String}
    doselection(options,Set(selector))
end

#selection for Strings
function doselection(options::Set{String},selector::String)::Vector{String}
    doselection(options,Set([selector]))
end

#do selection for Regex
function doselection(options::Set{String},selector::Regex)::Vector{String}
    selectedset=filter(options) do o
        !isnothing(match(selector,o))
    end |> Set
    doselection(options,selectedset)
end

#selection for Colon
function doselection(options::Set{String},::Colon)::Vector{String}
    collect(options)
end

"""
```julia
DeltaCT(data, hk; ntc="NTC", ampthreshold=40)
```

Perform ΔCT analysis on a `QPCRDataset` by normalizing the results for each sample 
against the geometric mean of the targets listed in `hk` (which should be a `Vector{String}`).

This function will throw and error if any CT values for the sample passed as `ntc` are greater
than `ampthreshold`.

# Indexing
The resulting object can be indexed using the form `result[targets]` or `result[targets,samples]`
where `targets` and `samples` can be a `String`, a `Vector{String}`, a `Regex` or `:`.
"""
struct DeltaCT <: DeltaResult
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
                error("Amplification detected in no template controls for the following targets: " * join(Set(amplified_ntc.target)," ,"))
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

    #one more inner constructor to allow for indexing of DeltaCT objects
    function DeltaCT(dct::DeltaCT,selectedtargets::Vector{String},selectedsamples::Vector{String})
        st=Set(selectedtargets)
        ss=Set(selectedsamples)
        #our hkgenes always get to come along for the ride
        st=union(st,Set(dct.hkgenes))
        newdata=selectdata(dct.data,st,ss)
        new(newdata,dct.hkgenes)
    end
end

"""
```julia
DDCT(dct, refsample)
```

Perform  ΔΔCT (comparitive CT) analysis on a `DeltaCT` object by renormalizing the ΔCT results
against the sample `refsample`

# Indexing
The resulting object can be indexed using the form `result[targets]` or `result[targets,samples]`
where `targets` and `samples` can be a `String`, a `Vector{String}`, a `Regex` or `:`.
"""
struct DDCT <: DeltaResult
    data::AbstractDataFrame
    dct::DeltaCT
    refsample::AbstractString
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

    #one more inner constructor to allow for indexing of DDCT objects
    function DDCT(ddct::DDCT,selectedtargets::Vector{String},selectedsamples::Vector{String})
        st=Set(selectedtargets)
        ss=Set(selectedsamples)
        #our hkgenes always get to come along for the ride
        st=union(st,Set(ddct.dct.hkgenes))
        #our refsample also are always included
        ss=union(ss,Set([ddct.refsample]))

        #also do the selection on our dct object
        new_dct=DeltaCT(ddct.dct,st |> collect,ss |> collect)

        #now select our dataframe
        newdata=selectdata(ddct.data,st |> collect,ss |> collect)
        new(newdata,new_dct,ddct.refsample)
    end
    
end

"""
```julia
writeresult(path,result)
result |> writeresult(path)
```
Write an analysis result to `path` as a CSV.
"""
function writeresult end

function writeresult(filename::AbstractString,result::DeltaResult)
    result.data |> CSV.write(filename)
end

function writeresult(filename::AbstractString)
    function(result::DeltaResult)
        writeresult(filename,result)
    end
end

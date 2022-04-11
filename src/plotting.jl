@recipe function plotdeltact(dct::DeltaCT)
    #remove the housekeeping genes, not interesting
    nohk=dct.data[map(dct.data.target) do t
                      !(t in dct.hkgenes)
                  end, :] |> copy
    #make a plot for each sample
    targetvec=Vector{Vector{String}}()
    dctfoldvec=Vector{Vector{Real}}()
    sampvec=String[]
    for s in Set(nohk.sample)
        push!(sampvec,s)
        this_samp=nohk[nohk.sample .== s,:]
        push!(targetvec,this_samp.target)
        push!(dctfoldvec,this_samp.dct_fold)
    end

    seriestype := :bar
    title := reshape(sampvec,1,:)
    xlabel := "Target"
    ylabel := "DCT Fold"
    legend --> false
    layout --> (length(sampvec),1)
    size --> (100*length(targetvec[1])+200,200*length(sampvec))
    left_margin --> 10mm
    (hcat(targetvec...),hcat(dctfoldvec...))
end

@recipe function plotddct(ddct::DDCT)
    #remove the housekeeping genes, not interesting
    nohk=ddct.data[map(ddct.data.target) do t
                      !(t in ddct.dct.hkgenes)
                  end, :] |> copy

    #also remove the refrence sample
    noref=nohk[map(nohk.sample) do s
                   s != ddct.refsample
               end, :] |> copy
    
    #make a plot for each target
    targetvec=Vector{String}()
    ddctfoldvec=Vector{Vector{Real}}()
    samplevec=Vector{Vector{String}}()
    for t in Set(noref.target)
        push!(targetvec,t)
        this_target=noref[noref.target .== t,:]
        push!(samplevec,this_target.sample)
        push!(ddctfoldvec,this_target.ddct_fold)
    end

    seriestype := :bar
    title := reshape(targetvec,1,:)
    xlabel := "Sample"
    ylabel := "DDCT Fold"
    legend --> false
    layout --> (length(targetvec),1)
    size --> (100*length(samplevec[1])+200,200*length(targetvec))
    left_margin --> 10mm
    (hcat(samplevec...),hcat(ddctfoldvec...))
end

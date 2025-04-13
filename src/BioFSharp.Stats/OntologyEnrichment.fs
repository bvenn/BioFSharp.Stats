namespace BioFSharp.Stats

open System
open FSharpAux

module OntologyEnrichment =

    /// Represents an item in an ontology set
    type OntologyItem<'a> = {
        Id               : string
        OntologyTerm     : string
        GroupIndex       : int
        Item             : 'a
        }


    /// Creates an item in an ontology set
    let createOntologyItem id ontologyTerm groupIndex item =
        {Id = id; OntologyTerm = ontologyTerm; GroupIndex = groupIndex; Item = item}

    /// Represents a gene set enrichment result 
    type GseaResult<'a> = {
        ///Ontology term e.g. MapMan term, GO term ...
        OntologyTerm     : string
        ///Sequence of single items associated to the ontology term 
        ItemsInBin       : seq<OntologyItem<'a>>
        ///Number of significantly altered items in 'OntologyTerm' bin
        NumberOfDEsInBin : int
        ///Number of items in 'OntologyTerm' bin
        NumberInBin      : int
        ///Number of significantly altered items within the total data set
        TotalNumberOfDE  : int
        ///Number of all items (expanded)
        TotalUniverse    : int
        ///p value as calculated by hypergeometric test
        PValue           : float
        } with
        static member Create(ontologyTerm, itemsInBin, numberOfDEsInBin, numberInBin, totalNumberOfDE, totalUnivers, pValue, ?FDR: float, ?QVal: float, ?Pi0: float) =
            {OntologyTerm = ontologyTerm; ItemsInBin = itemsInBin; NumberOfDEsInBin = numberOfDEsInBin; 
                NumberInBin = numberInBin; TotalNumberOfDE = totalNumberOfDE; TotalUniverse = totalUnivers; PValue = pValue;}

    /// Represents a gene set enrichment result with extended statistics 
    type GseaResultFdrExtended<'a> = {
        ///Ontology term e.g. MapMan term, GO term ...
        OntologyTerm     : string
        ///Sequence of single items associated to the ontology term 
        ItemsInBin       : seq<OntologyItem<'a>>
        ///Number of significantly altered items in 'OntologyTerm' bin
        NumberOfDEsInBin : int
        ///Number of items in 'OntologyTerm' bin
        NumberInBin      : int
        ///Number of significantly altered items within the total data set
        TotalNumberOfDE  : int
        ///Number of all items (expanded)
        TotalUniverse    : int
        ///p value as calculated by hypergeometric test
        PValue           : float
        ///false discovery rate (FDR) using Benjamini-Hochberg method
        FDR_BH           : float
        ///q value (FDR) as calculated by Storey method
        QValue           : float
        ///pi0 (proportion of null hypotheses) as calculated by Storey method
        Pi0              : float
        } with
        static member Create(ontologyTerm, itemsInBin, numberOfDEsInBin, numberInBin, totalNumberOfDE, totalUnivers, pValue, fdr: float, qVal: float, pi0: float) =
            {OntologyTerm = ontologyTerm; ItemsInBin = itemsInBin; NumberOfDEsInBin = numberOfDEsInBin; 
                NumberInBin = numberInBin; TotalNumberOfDE = totalNumberOfDE; TotalUniverse = totalUnivers; PValue = pValue; FDR_BH = fdr; QValue = qVal; Pi0 = pi0}

    ///Splits an OntologyEntry with seperator concatenated TermIds
    let splitMultipleAnnotationsBy (separator:char) (item:OntologyItem<'A>) =
        let annotations = item.OntologyTerm.Split(separator)
        annotations
        |> Seq.map (fun ot -> {item with OntologyTerm = ot})

    /// Splits MapMan OntologyEntries with seperator concatenated TermIds
    /// Attention: Also parses string to int to get rid of 0 - terms
    let splitMapManOntologyItemsBy (separator:char) (data:seq<OntologyItem<'a>>) =
        let splitTerm (termId:string) (separator:char) =
            termId.Split(separator) 
            |> Array.map (fun sTerm -> 
                let splited = sTerm.Split('.')
                let toInt = splited |> Seq.map (fun v -> Int32.Parse(v).ToString())                                                                            
                toInt  |> String.concat "." 
                         )
        data
        |> Seq.collect (fun oi -> 
            splitTerm oi.OntologyTerm separator
            |> Seq.map (fun sTerm -> createOntologyItem oi.Id sTerm oi.GroupIndex oi.Item)
                        )


    /// Extends leaf OntologyEntries to their full tree
    let expandOntologyTree (data:seq<OntologyItem<'a>>) =
        data
        |> Seq.collect (fun oi -> 
            let expandenTermIds = oi.OntologyTerm.Split('.') |> Array.scanReduce (fun acc elem -> acc + "." + elem)
            expandenTermIds |> Seq.map (fun sTerm -> createOntologyItem oi.Id sTerm oi.GroupIndex oi.Item) 
                       )

    /// <summary>Calculates p value based on hypergeometric distribution (pValue <= k)</summary>
    /// <remarks>the hypergeometric distribution is a discrete probability distribution that describes the probability of 
    ///   k successes in
    ///   n draws from a finite 
    ///   N population of size containing
    ///   K successes without replacement (successes states)</remarks>
    /// <param name="numberOfDEsInBin">(k) number of significantly altered entities in the respective 'OntologyTerm' bin</param>
    /// <param name="numberInBin">(n) number of entities within the respective 'OntologyTerm' bin</param>
    /// <param name="totalUnivers">(N) total number of entities in the experiment</param>
    /// <param name="totalNumberOfDE">(K) total number of significantly altered entities in the dataset</param>
    /// <param name="splitPvalueThreshold">threshold until mid p values are calculated (default=5)</param>
    /// <returns>A pvalue determined by hypergeometric test. If the number of entities within a bin is below the splitPValue
    /// threshold a mid p value is calculated that is lower than its default</returns>
    let calcHyperGeoPvalue numberOfDEsInBin numberInBin totalUnivers totalNumberOfDE (splitPvalueThreshold:int) =
        let hp = FSharp.Stats.Distributions.Discrete.Hypergeometric.Init totalUnivers totalNumberOfDE numberInBin            
        if numberInBin > splitPvalueThreshold then                                
            // Calculate normal pValue
            1. -  hp.CDF (float (numberOfDEsInBin - 1)) 
        else
            // Calculate split pValue
            0.5 * ((1. -  hp.CDF(float(numberOfDEsInBin - 1)) ) + ( (1. -  hp.CDF(float(numberOfDEsInBin))) ) )

    
    /// <summary>Calculates functional term enrichment</summary>
    /// <remarks>http://bioinformatics.oxfordjournals.org/cgi/content/abstract/23/4/401</remarks>
    /// <param name="deGroupIndex">defines the index of the differential expression group</param>
    /// <param name="splitPvalueThreshold">threshold until mid p values are calculated (default=5)</param>
    /// <param name="minNumberInTerm">describes how many items are required in a bin to be taken into account for multiple testing correction (default = 2)</param>
    /// <param name="data">sequence of ontology items</param>
    /// <returns>sequence of GseaResult</returns>
    /// <example>
    /// <code>
    /// let data = [|
    ///     createOntologyItem "id1" "photosynthesis.lightreaction" 0 "item1"
    ///     createOntologyItem "id2" "protein.degradation" 0 "item2"
    ///    createOntologyItem "id3" "photosynthesis.lightreaction.LHCI" 1 "item3" // significantly different gene/protein/metabolite
    /// calcOverEnrichment 1 (Some 5) (Some 2) data
    /// </code>
    /// </example>
    let calcOverEnrichment (deGroupIndex:int) (splitPvalueThreshold:option<int>) (minNumberInTerm:option<int>) (data:seq<OntologyItem<'a>>) =
        let _splitPvalueThreshold   = defaultArg splitPvalueThreshold 0
        let _minNumberInBin        = defaultArg minNumberInTerm 0
        
        // Distinct by term and gene name
        // Has to be done by an ouside function
        //let distinctData    = data |> Seq.distinctBy (fun o -> o.displayID)                
        let gData           = data |> Seq.groupBy ( fun o -> o.OntologyTerm)
        // reduce to terms at least annotated with 2 items
        let fData = gData |> Seq.filter ( fun (key:string,values:seq<OntologyItem<'a>>) -> Seq.length(values) >= _minNumberInBin)
        let groupCount = fData |> Seq.collect (fun (key:string,values:seq<OntologyItem<'a>>) -> values ) |> Seq.countBy (fun o -> o.GroupIndex)
        
        let totalUniverse    = groupCount |> Seq.fold (fun  (acc:int) (index:int,count:int) -> acc + count) 0
        let totalNumberOfDE = 
            let tmp = groupCount |> Seq.tryFind (fun (key,v) -> key = deGroupIndex)
            if tmp.IsNone then 
                raise (System.ArgumentException("DE group index does not exists in ontology entry"))
            else
                snd(tmp.Value)
        
        // returns (DE count, all count)
        let countDE (subSet:seq<OntologyItem<'a>>) =             
            let countMap = 
                subSet 
                |> Seq.countBy (fun (oi) -> oi.GroupIndex = deGroupIndex)
                |> Map.ofSeq
            (countMap.TryFindDefault 0 true,(countMap.TryFindDefault 0 true) + (countMap.TryFindDefault 0 false))
        
        let results = 
            fData
            |> Seq.map (fun (oTerm,values) -> 
                let numberOfDEsInBin,numberInBin = countDE values
                let pValue = calcHyperGeoPvalue numberOfDEsInBin numberInBin totalUniverse totalNumberOfDE _splitPvalueThreshold
                GseaResult<'a>.Create(oTerm, values, numberOfDEsInBin, numberInBin, totalNumberOfDE, totalUniverse, pValue))
        
        results

    [<Obsolete("Use OntologyEnrichment.calcOverEnrichment instead")>]
    let CalcOverEnrichment (deGroupIndex:int) (splitPvalueThreshold:option<int>) (minNumberInTerm:option<int>) (data:seq<OntologyItem<'a>>)= 
        calcOverEnrichment deGroupIndex splitPvalueThreshold minNumberInTerm data 

    [<Obsolete("Use OntologyEnrichment.calcOverEnrichment instead")>]
    let CalcSimpleOverEnrichment (deGroupIndex:int) (splitPvalueThreshold:option<int>) (minNumberInTerm:option<int>) (data:seq<OntologyItem<'a>>)= 
        calcOverEnrichment deGroupIndex splitPvalueThreshold (Some 0) data 

    /// <summary>Calculates functional term enrichment with additional multiple testing correction</summary>
    /// <remarks>http://bioinformatics.oxfordjournals.org/cgi/content/abstract/23/4/401, storey q value, benjamini hochberg FDR</remarks>
    /// <param name="deGroupIndex">defines the index of the differential expression group</param>
    /// <param name="splitPvalueThreshold">threshold until mid p values are calculated (default=5)</param>
    /// <param name="minNumberInTerm">describes how many items are required in a bin to be taken into account for multiple testing correction (default = 2)</param>
    /// <param name="data">sequence of ontology items</param>
    /// <returns>sequence of GseaResultFdrExtended</returns>
    /// <example>
    /// <code>
    /// let data = [|
    ///     createOntologyItem "id1" "photosynthesis.lightreaction" 0 "item1"
    ///     createOntologyItem "id2" "protein.degradation" 0 "item2"
    ///    createOntologyItem "id3" "photosynthesis.lightreaction.LHCI" 1 "item3" // significantly different gene/protein/metabolite
    /// calcOverEnrichmentIncludeFDR 1 (Some 5) (Some 2) data
    /// </code>
    /// </example>
    let calcOverEnrichmentIncludeFDR (deGroupIndex:int) (splitPvalueThreshold:option<int>) (minNumberInTerm:option<int>) (data:seq<OntologyItem<'a>>) =
        let results = 
            calcOverEnrichment deGroupIndex splitPvalueThreshold minNumberInTerm data 

        let pvalues = 
            results 
            |> Seq.map _.PValue 
            |> Array.ofSeq
    
        // fdr calculation 
        let fdr_bh = 
            FSharp.Stats.Testing.MultipleTesting.benjaminiHochbergFDR pvalues 
            |> Array.ofSeq
        
        // q value calculation
        let pi0 = FSharp.Stats.Testing.MultipleTesting.Qvalues.pi0Bootstrap pvalues
        let qval = FSharp.Stats.Testing.MultipleTesting.Qvalues.ofPValues pi0 pvalues
        results
        |> Seq.mapi (fun i res -> 
            let fdr = fdr_bh.[i]
            let qval = qval.[i]
            let pi0 = pi0
            GseaResultFdrExtended<'a>.Create(
                res.OntologyTerm, 
                res.ItemsInBin, 
                res.NumberOfDEsInBin, 
                res.NumberInBin, 
                res.TotalNumberOfDE, 
                res.TotalUniverse, 
                res.PValue, 
                fdr, qval, pi0)
            
        )
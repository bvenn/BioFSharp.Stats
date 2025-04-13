namespace BioFSharp.Stats.Tests

open Xunit
open BioFSharp
open BioFSharp.Stats
open OntologyEnrichment

module OntologyEnrichmentTests =

    let data =  [|
        createOntologyItem "id01" "PS.lightreaction.LHCII;protein.degradation.ubiquitin" 1 "description01"
        createOntologyItem "id02" "PS.lightreaction.LHCII" 0 "description02"
        createOntologyItem "id03" "PS.lightreaction.LHCII" 1 "description03"
        createOntologyItem "id04" "PS.lightreaction.LHCII" 0 "description04"
        createOntologyItem "id05" "PS.lightreaction.LHCII" 0 "description05"
        createOntologyItem "id06" "PS.lightreaction.LHCII" 1 "description06"
        createOntologyItem "id07" "PS.lightreaction.LHCII" 0 "description07"
        createOntologyItem "id08" "PS.lightreaction.LHCII" 1 "description08"
        createOntologyItem "id09" "PS.lightreaction.LHCII" 0 "description09"
        createOntologyItem "id10" "PS.lightreaction.LHCII" 1 "description10"
        createOntologyItem "id11" "PS.lightreaction.LHCII" 0 "description11" //36
        createOntologyItem "id12" "PS.lightreaction" 0 "description12"
        createOntologyItem "id13" "PS.lightreaction" 0 "description13"
        createOntologyItem "id14" "PS.lightreaction" 0 "description14"
        createOntologyItem "id15" "PS.lightreaction" 1 "description15" //8
        createOntologyItem "id16" "protein.degradation.ubiquitin" 0 "description16"
        createOntologyItem "id17" "protein.degradation.ubiquitin" 1 "description17"
        createOntologyItem "id18" "protein.degradation.ubiquitin" 0 "description18"
        createOntologyItem "id19" "protein.degradation.ubiquitin" 0 "description19"
        createOntologyItem "id20" "protein.degradation.ubiquitin" 1 "description20"
        createOntologyItem "id21" "protein.degradation.ubiquitin" 0 "description21" //18
        createOntologyItem "id22" "protein.degradation.ubiquitin.e1" 0 "description22"
        createOntologyItem "id23" "protein.degradation.ubiquitin.e1" 0 "description23"
        createOntologyItem "id24" "protein.degradation.ubiquitin.e1" 0 "description24" //12
        createOntologyItem "id25" "protein.synthesis.initiation" 1 "description25"
        createOntologyItem "id26" "protein.synthesis.initiation" 0 "description26"
        createOntologyItem "id27" "protein.synthesis.initiation" 1 "description27"
        createOntologyItem "id28" "protein.synthesis.initiation" 1 "description28"
        createOntologyItem "id29" "protein.synthesis.initiation" 1 "description29"
        createOntologyItem "id30" "protein.synthesis" 0 "description30" //17
        createOntologyItem "id31" "singletest" 1 "description31" //1
    |]

    let dataSplit = 
        data
        |> Seq.collect (splitMultipleAnnotationsBy ';')

    let dataExtended = 
        expandOntologyTree dataSplit

    let enrichmentResult = 
        calcOverEnrichmentIncludeFDR 1 (Some 0) (Some 0) dataExtended

    [<Fact>]
    let ``splitOntologyItems`` () =
        let actual = dataSplit |> Seq.length
        let expected = data.Length + 1 //first one is split into 2, all others remain the same
        Assert.Equal<int>(expected, actual)

    [<Fact>]
    let ``expandOntologyTree`` () =
        let actual = dataExtended |> Seq.length
        let expected = 92 // 36 + 8 + 18 + 12 + 18
        Assert.Equal<int>(expected, actual)

    [<Fact>]
    let ``ontologyEnrichment`` () =
        // per bin one result
        let numberOfBins_actual = enrichmentResult |> Seq.length
        let numberOfBins_expected = 10
        Assert.Equal<int>(numberOfBins_expected, numberOfBins_actual)
        let bin1 = enrichmentResult |> Seq.find (fun x -> x.OntologyTerm = "PS.lightreaction.LHCII")
        Assert.Equal<int>(11, bin1.NumberInBin)
        Assert.Equal<int>(5, bin1.NumberOfDEsInBin)
        Assert.Equal<int>(39, bin1.TotalNumberOfDE)
        Assert.Equal<int>(92, bin1.TotalUniverse)
        Assert.Equal<float>(0.5367813329, System.Math.Round(bin1.PValue, 10)) //https://systems.crump.ucla.edu/hypergeometric/index.php
        let testbin = 
            enrichmentResult |> Seq.find (fun x -> x.OntologyTerm = "singletest")
        // test for a singleton bin of which valid p values exist
        Assert.Equal<float>(0.4239130435, System.Math.Round(testbin.PValue, 10))
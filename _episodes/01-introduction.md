---
title: "General Physics Objects and POET"
teaching: 20
exercises: 40
questions:
- What do we call physics objects in CMS?
- How are physics objects reconstructed?
- What is POET?
- How do I run POET?
objectives:
- Learn about the different physics objects in CMS and get briefed on their reconstruction
- Learn what the POET software is and how to run it
- Learn the structure of the POET configurator and how add new physics objects.
keypoints:
- Physics objects are the final abstraction in the detector that can be associated to physical entities like particles.
- POET is a collection of CMSSW EDAnalyzers meant to teach how to access physics objects information.
---

## Overview

The **CMS experment** is a giant detector that acts like a camera that "photographs" particle collisions, allowing us to interpret their nature.

Certainly, we cannot directly observe all the particles created in the collisions because some of them decay very quickly or simply do not interact with our detector. However, *we can infer their presence*. If they decay to other stable particles and interact with the apparatus, they leave signals in the CMS subdetectors. These signals are used to reconstruct the decay products or infer their presence; we call these, **physics objects**. 

Physics objects are built with the information collected by the sensors of the CMS detector.  Take a look at the CMS experiment in the image below. This is a good recent [video](https://youtu.be/6gfvGTCWXaw) that you can watch later to get a quick feeling of how CMS looks now in Run 3.

<img src="../fig/cmsdetector.png" width="800">

Physics objects could be

*   muons
*   electrons
*   jets
*   photons
*   taus
*   missing transverse momentum



For the current releases of open data, we store them in ROOT files following the EDM data model in the so-called **miniAOD** format.  One needs the **CMSSW** software to read and extract information (of the physics objects) from these files.

In the [CERN Open Portal](http://opendata.cern.ch) (CODP) site one can find a more detailed description of these physical objects and a list of them corresponding to [2010](http://opendata.cern.ch/docs/cms-physics-objects-2010), [2011/2012](http://opendata.cern.ch/docs/cms-physics-objects-2011) from Run 1, and [2015](http://opendata.cern.ch/docs/cms-physics-objects-2015) (Run 2) releases of open data.

In this workshop we will focus on working with open data from the **latest 2015 release** from Run 2. 


## Physics Objects reconstruction

Physics objects are mainly **reconstructed** using methods like clustering and linking parts of a CMS subsystem to parts of other CMS subsystems.  For instance, electromagnetic objects (electrons or photons) are reconstructed by the linking of tracks with ECAL energy deposits.

<img src="../fig/cmsslice.png" width="1216">


These actions are essential parts of the so-called **Particle Flow** (PF) algorithm.


<img src="../fig/pflow.png" width="600">

The **particle-flow (PF) algorithm** aims at reconstructing and identifying all stable particles in the event, i.e., electrons, muons, photons, charged hadrons and neutral hadrons, with a thorough combination of all CMS sub-detectors towards an optimal determination of their direction, energy and type. This **list of individual particles is then used**, as if it came from a Monte-Carlo event generator, to build jets (from which the quark and gluon energies and directions are inferred), to determine the missing transverse energy (which gives an estimate of the direction and energy of the neutrinos and other invisible particles), to reconstruct and identify taus from their decay products and more.

The PF algorithm and all its reconstruction tributaries are written in C++ and are part of the **CMSSW** software. The PF algorithm will be covered in detail during the [Advanced Tools lesson](https://cms-opendata-workshop.github.io/workshop2023-advobjects/) of the live workshop.

## The Physics Object Extractor Tool (POET)

The 2015 POET repository instructions read *"... POET repository contains packages that provide instructions and examples on how to extract physics (objects) information from Run 2 (MiniAOD format) CMS open/legacy data and methods or tools needed for processing them"*

POET was thought, originally, as tool to **learn** how to use CMSSW to access objects in an easier fashion.  The code written by CMS experienced users can be extremely convoluted, but it all follows the same logic.  POET tries to show how to do it in a pedagogical way by separating the processing of different objects.  POET is just a **collection of CMSSW EDAnalyzers**, which you could learn about by follow a [lesson from the 2023 Open Data Workshop](https://cms-opendata-workshop.github.io/workshop2023-lesson-cmssw/). It is also similar in structure and content to ["NanoAOD"](https://github.com/cms-opendata-analyses/AOD2NanoAODOutreachTool), a CMS data format developed during Run 2 that is even smaller than MiniAOD. 

For this short workshop, our goal is to give you a flavor of how physics object information is stored, and how the POET can be run to produce ROOT files in a simple format for analysts to process. **Most users of Open Data will not need to edit or change the configuration of POET**: for any analysis that can be built using POET-like (or ``NanoAOD-like") ROOT files, POET could be run as a first step in a scaled-up analysis workflow. 

### Downloading POET

As a useful reference, let us start with learning how to run with POET.

First, fire up your [`CMSSW` Docker container](https://cms-opendata-workshop.github.io/workshopwhepp-lesson-docker/03-docker-for-cms-opendata/index.html#download-the-docker-image-for-cmssw-open-data-and-start-a-container) already used during the pre-exercises. If you have the Mac M1/M2 chip, this code could be downloaded and investigated (but not run) in the ROOT container.

~~~
docker start -i my_od #use the name you gave to yours
~~~ 
{: .language-bash}

Once inside your container, make sure you are at the `src` directory level of your `CMSSW` area
~~~
pwd
~~~ 
{: .language-bash}

~~~
/code/CMSSW_7_6_7/src
~~~ 
{: .output}

Let's clone the POET brach that we will use for this lesson:

~~~
git clone -b 2015MiniAOD https://github.com/cms-opendata-analyses/PhysObjectExtractorTool.git
~~~ 
{: .language-bash}

### Exploring the POET setup (Optional)

Now, `cd` to the package of interest:

~~~
cd PhysObjectExtractorTool/PhysObjectExtractor/
~~~ 
{: .language-bash}

Explore this directory with a `ls` command.  
Explore the `src` directory, you will find several *EDAnalyzers*, essentially one for each object of interest.

~~~
ls src
~~~ 
{: .language-bash}

~~~
ElectronAnalyzer.cc  GenParticleAnalyzer.cc  MetAnalyzer.cc   PhotonAnalyzer.cc     TauAnalyzer.cc          TriggerAnalyzer.cc
FatjetAnalyzer.cc    JetAnalyzer.cc          MuonAnalyzer.cc  SimpleEleMuFilter.cc  TriggObjectAnalyzer.cc  VertexAnalyzer.cc
~~~ 
{: .output}

Explore the `python` directory, and you will find some configuration files. 
There is one master configuration file to configure anything that we want to do with POET: `poet_cfg.py`. 
To begin with, let's study a short example, `poet_demo_cfg.py`. Open `poet_demo_cfg.py` using a text editor, and this is what you will see:

~~~
import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes
import sys

#---- sys.argv takes the parameters given as input cmsRun PhysObjectExtractor/python/poet_cfg.py <isData (default=False)>
#----  e.g: cmsRun PhysObjectExtractor/python/poet_cfg.py True
#---- NB the first two parameters are always "cmsRun" and the config file name
#---- Work with data (if False, assumed MC simulations)
#---- This needs to be in agreement with the input files/datasets below.
if len(sys.argv) > 2:
    isData = eval(sys.argv[2])
else:
    isData = False
isMC = True
if isData: isMC = False

process = cms.Process("POET")

#---- Configure the framework messaging system
#---- https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMessageLogger
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = "WARNING"
process.MessageLogger.categories.append("POET")
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit=cms.untracked.int32(-1))
process.options = cms.untracked.PSet(wantSummary=cms.untracked.bool(True))

#---- Select the maximum number of events to process (if -1, run over all events)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

#---------------------------------------------------- 
#---- More configuration skipped here in the webpage
#----     |
#----     |
#----     v
#---------------------------------------------------- 

#---- Define the test source files to be read using the xrootd protocol (root://), or local files (file:)
process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(
        #'root://eospublic.cern.ch//eos/opendata/cms/mc/RunIIFall15MiniAODv2/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/00000/0005EA25-8CB8-E511-A910-00266CF85DA0.root'   
        'root://eospublic.cern.ch//eos/opendata/cms/mc/RunIIFall15MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext3-v1/00000/02837459-03C2-E511-8EA2-002590A887AC.root'
        )
)
if isData:
    process.source.fileNames = cms.untracked.vstring(
        #'root://eospublic.cern.ch//eos/opendata/cms/Run2015D/DoubleMuon/MINIAOD/16Dec2015-v1/10000/000913F7-E9A7-E511-A286-003048FFD79C.root'
 'root://eospublic.cern.ch//eos/opendata/cms/Run2015D/SingleElectron/MINIAOD/08Jun2016-v1/10000/001A703B-B52E-E611-BA13-0025905A60B6.root'
        )


#----- Configure POET analyzers -----#
process.myelectrons = cms.EDAnalyzer('ElectronAnalyzer', 
                                     electrons = cms.InputTag("slimmedElectrons"), 
                                     vertices=cms.InputTag("offlineSlimmedPrimaryVertices"))

process.mymuons = cms.EDAnalyzer('MuonAnalyzer', 
                                 muons = cms.InputTag("slimmedMuons"), 
                                 vertices=cms.InputTag("offlineSlimmedPrimaryVertices"))




#----- RUN THE JOB! -----#
process.TFileService = cms.Service("TFileService", fileName=cms.string("myoutput.root"))

if isData:
 process.p = cms.Path(process.myelectrons)
else:
 process.p = cms.Path(process.myelectrons)
~~~ 
{: .language-python}


After important several key CMS software packages, the configuration file can interpret an argument from the user:

~~~
#---- sys.argv takes the parameters given as input cmsRun PhysObjectExtractor/python/poet_cfg.py <isData (default=False)>
#----  e.g: cmsRun PhysObjectExtractor/python/poet_cfg.py True
#---- NB the first two parameters are always "cmsRun" and the config file name
#---- Work with data (if False, assumed MC simulations)
#---- This needs to be in agreement with the input files/datasets below.
if len(sys.argv) > 2:
    isData = eval(sys.argv[2])
else:
    isData = False
isMC = True
if isData: isMC = False

process = cms.Process("POET")
~~~ 
{: .language-python}

This basically tells you that you have the option of choosing between running over collisions data (or just data) or over montecarlo (MC) simulations.  If you don't give any argument, the `isMC` will remain true, whereas if you pass a `True` boolean, then the `isData` switch will become active.

This will have a direct impact on the kind of input file the `process.source` module will read, and whether (for collision data) the data quality filter will be applied. For convinience, we have put the .json file for the quality filter in the `data` directory of the POET repository. More information can be found [here](https://opendata.cern.ch/record/14210).

~~~
#---- Define the test source files to be read using the xrootd protocol (root://), or local files (file:)
process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(
        #'root://eospublic.cern.ch//eos/opendata/cms/mc/RunIIFall15MiniAODv2/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/00000/0005EA25-8CB8-E511-A910-00266CF85DA0.root'   
        'root://eospublic.cern.ch//eos/opendata/cms/mc/RunIIFall15MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext3-v1/00000/02837459-03C2-E511-8EA2-002590A887AC.root'
        )
)
if isData:
    process.source.fileNames = cms.untracked.vstring(
        #'root://eospublic.cern.ch//eos/opendata/cms/Run2015D/DoubleMuon/MINIAOD/16Dec2015-v1/10000/000913F7-E9A7-E511-A286-003048FFD79C.root'
 'root://eospublic.cern.ch//eos/opendata/cms/Run2015D/SingleElectron/MINIAOD/08Jun2016-v1/10000/001A703B-B52E-E611-BA13-0025905A60B6.root'
        )

    #---- Apply the data quality JSON file filter. This example is for 2015 data
    #---- It needs to be done after the process.source definition
    #---- Make sure the location of the file agrees with your setup
    goodJSON = "data/Cert_13TeV_16Dec2015ReReco_Collisions15_25ns_JSON_v2.txt"
    myLumis = LumiList.LumiList(filename=goodJSON).getCMSSWString().split(",")
    process.source.lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange())
    process.source.lumisToProcess.extend(myLumis) 

~~~ 
{: .language-python}

Let's explore one of these input files. These are the data we release under the open data initiative.  You notice from their name that their format is `MINIAOD`.  A given dataset contains many of these files.  Choose any you like and check its content, for instance:

~~~
edmDumpEventContent root://eospublic.cern.ch//eos/opendata/cms/mc/RunIIFall15MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext3-v1/00000/02837459-03C2-E511-8EA2-002590A887AC.root
~~~ 
{: .language-bash}

~~~
Type                                  Module                      Label             Process   
----------------------------------------------------------------------------------------------
LHEEventProduct                       "externalLHEProducer"       ""                "LHE"     
GenEventInfoProduct                   "generator"                 ""                "SIM"     
edm::TriggerResults                   "TriggerResults"            ""                "SIM"     
edm::TriggerResults                   "TriggerResults"            ""                "HLT"     
HcalNoiseSummary                      "hcalnoise"                 ""                "RECO"    
L1GlobalTriggerReadoutRecord          "gtDigis"                   ""                "RECO"    
double                                "fixedGridRhoAll"           ""                "RECO"    
double                                "fixedGridRhoFastjetAll"    ""                "RECO"    
double                                "fixedGridRhoFastjetAllCalo"   ""                "RECO"    
double                                "fixedGridRhoFastjetCentral"   ""                "RECO"    
double                                "fixedGridRhoFastjetCentralCalo"   ""                "RECO"    
double                                "fixedGridRhoFastjetCentralChargedPileUp"   ""                "RECO"    
double                                "fixedGridRhoFastjetCentralNeutral"   ""                "RECO"    
edm::TriggerResults                   "TriggerResults"            ""                "RECO"    
reco::BeamHaloSummary                 "BeamHaloSummary"           ""                "RECO"    
reco::BeamSpot                        "offlineBeamSpot"           ""                "RECO"    
vector<l1extra::L1EmParticle>         "l1extraParticles"          "Isolated"        "RECO"    
vector<l1extra::L1EmParticle>         "l1extraParticles"          "NonIsolated"     "RECO"    
vector<l1extra::L1EtMissParticle>     "l1extraParticles"          "MET"             "RECO"    
vector<l1extra::L1EtMissParticle>     "l1extraParticles"          "MHT"             "RECO"    
vector<l1extra::L1HFRings>            "l1extraParticles"          ""                "RECO"    
vector<l1extra::L1JetParticle>        "l1extraParticles"          "Central"         "RECO"    
vector<l1extra::L1JetParticle>        "l1extraParticles"          "Forward"         "RECO"    
vector<l1extra::L1JetParticle>        "l1extraParticles"          "IsoTau"          "RECO"    
vector<l1extra::L1JetParticle>        "l1extraParticles"          "Tau"             "RECO"    
vector<l1extra::L1MuonParticle>       "l1extraParticles"          ""                "RECO"    
edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> >    "reducedEgamma"             "reducedEBRecHits"   "PAT"     
edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> >    "reducedEgamma"             "reducedEERecHits"   "PAT"     
edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> >    "reducedEgamma"             "reducedESRecHits"   "PAT"     
edm::TriggerResults                   "TriggerResults"            ""                "PAT"     
edm::ValueMap<float>                  "offlineSlimmedPrimaryVertices"   ""                "PAT"     
pat::PackedTriggerPrescales           "patTrigger"                ""                "PAT"     
pat::PackedTriggerPrescales           "patTrigger"                "l1max"           "PAT"     
pat::PackedTriggerPrescales           "patTrigger"                "l1min"           "PAT"     
vector<PileupSummaryInfo>             "slimmedAddPileupInfo"      ""                "PAT"     
vector<pat::Electron>                 "slimmedElectrons"          ""                "PAT"     
vector<pat::Jet>                      "slimmedJets"               ""                "PAT"     
vector<pat::Jet>                      "slimmedJetsAK8"            ""                "PAT"     
vector<pat::Jet>                      "slimmedJetsPuppi"          ""                "PAT"     
vector<pat::Jet>                      "slimmedJetsAK8PFCHSSoftDropPacked"   "SubJets"         "PAT"     
vector<pat::Jet>                      "slimmedJetsCMSTopTagCHSPacked"   "SubJets"         "PAT"     
vector<pat::MET>                      "slimmedMETs"               ""                "PAT"     
vector<pat::MET>                      "slimmedMETsPuppi"          ""                "PAT"     
vector<pat::Muon>                     "slimmedMuons"              ""                "PAT"     
vector<pat::PackedCandidate>          "lostTracks"                ""                "PAT"     
vector<pat::PackedCandidate>          "packedPFCandidates"        ""                "PAT"     
vector<pat::PackedGenParticle>        "packedGenParticles"        ""                "PAT"     
vector<pat::Photon>                   "slimmedPhotons"            ""                "PAT"     
vector<pat::Tau>                      "slimmedTaus"               ""                "PAT"     
vector<pat::TriggerObjectStandAlone>    "selectedPatTrigger"        ""                "PAT"     
vector<reco::CATopJetTagInfo>         "caTopTagInfosPAT"          ""                "PAT"     
vector<reco::CaloCluster>             "reducedEgamma"             "reducedEBEEClusters"   "PAT"     
vector<reco::CaloCluster>             "reducedEgamma"             "reducedESClusters"   "PAT"     
vector<reco::Conversion>              "reducedEgamma"             "reducedConversions"   "PAT"     
vector<reco::Conversion>              "reducedEgamma"             "reducedSingleLegConversions"   "PAT"     
vector<reco::GenJet>                  "slimmedGenJets"            ""                "PAT"     
vector<reco::GenJet>                  "slimmedGenJetsAK8"         ""                "PAT"     
vector<reco::GenParticle>             "prunedGenParticles"        ""                "PAT"     
vector<reco::GsfElectronCore>         "reducedEgamma"             "reducedGedGsfElectronCores"   "PAT"     
vector<reco::PhotonCore>              "reducedEgamma"             "reducedGedPhotonCores"   "PAT"     
vector<reco::SuperCluster>            "reducedEgamma"             "reducedSuperClusters"   "PAT"     
vector<reco::Vertex>                  "offlineSlimmedPrimaryVertices"   ""                "PAT"     
vector<reco::VertexCompositePtrCandidate>    "slimmedSecondaryVertices"   ""                "PAT"
~~~ 
{: .output}


This is essentially all the information you can get from MINIAOD files like this one. Put special attention to the anything related to *Electron(s)* or *Muon(s)*.  One thing you notice is that the acronym **PAT** is repeated several times.  PAT stands for Physics Analysis Toolkit and is a framework within CMSSW that is extensively used to refine the selection of physical objects in CMS.  The **RECO** variables are, however, lower level variables.  There are also some objects that will give you additional information, like the **HLT** *TriggerResults*, which we will cover separately.


Next in the `poet_demo_cfg.py` file, you will find a couple of modules that are evidently associated to EDAnalyzers that deal with electrons and muons:

~~~
#----- Configure POET analyzers -----#
process.myelectrons = cms.EDAnalyzer('ElectronAnalyzer', 
                                     electrons = cms.InputTag("slimmedElectrons"), 
                                     vertices=cms.InputTag("offlineSlimmedPrimaryVertices"))

process.mymuons = cms.EDAnalyzer('MuonAnalyzer', 
                                 muons = cms.InputTag("slimmedMuons"), 
                                 vertices=cms.InputTag("offlineSlimmedPrimaryVertices"))
~~~ 
{: .language-python}


The construction `cms.EDAnalyzer('ElectronAnalyzer'`, automatically tells you that there must be an EDAnalyzer called `ElectronAnalyzer` somewhere in this package.  Indeed, if you list your `src` directory, you will find, among other EDAnalyzers, an `ElectronAnalyzer.cc` file.  So, this module `process.myelectrons` configures that C++ code.   Something similar is true for the `MuonAnalyzer`.  Also notice that `InputTag`s used in these modules correspond to the collections seen above when we dump the structure of the input file.

Finally, there is a `TFileService` module.  This will cause CMSSW to create a `ROOT` file called `myoutput.root`, which will contain the information that a given EDAnalyzer returns.

~~~
process.TFileService = cms.Service("TFileService", fileName=cms.string("myoutput.root"))
~~~ 
{: .language-python}


> ### Enough talk, let's run POET! (Optional)
>
> Compile the code with `scram b` and then run POET twice with the `True` argument:
> 
> * First with `cmsRun python/poet_demo_cfg.py True` 
> * and then with just `cmsRun python/poet_demo_cfg.py` (without the *True* boolean, so we run over MC simulations)
>
> Open the `myoutput.root` file that gets produced using `root -l myoutput.root` and have a quick look with the `TBrowser`.  
>
> You will notice that only `electron` variables got stored.  Can you fix the `poet_demo_cfg.py` so we get the information from muons as well?
>
> > ## Solution
> >
> > Just add the `process.mymuons` module to the running paths at the end of the config file. Use a `+` sign to connect the modules in a sequence.  Run it again and check if you get some muon info (you don't need to recompile, that is the beauty of config files).
> > 
> {: .solution}
{: .challenge}

Congratulations, you have processed a CMS MiniAOD file and stored some physics objects!


### POET structure of EDAnalyzers (Optional)

We will briefly study the structure of an EDAnalyzer in POET, using `src/ElectronAnalyzer.cc` as an example.
The first thing to to notice is that all the EDAnalyzers in POET book a `ROOT` *TTree* and some variables in the class definition to store the physics objects information:

~~~
      TTree *mtree;
      int numelectron; //number of electrons in the event
      std::vector<float> electron_e;
      std::vector<float> electron_pt;
      std::vector<float> electron_px;
      std::vector<float> electron_py;
      std::vector<float> electron_pz;
      std::vector<float> electron_eta;
      std::vector<float> electron_phi;
~~~
{: .language-cpp}

In the constructor, a *TFileService* object is used to deal with the final `ROOT` output file and some `ROOT` braches are specified.  Of course, these are the ones which appear in the `myoutput.root` file.

~~~
   edm::Service<TFileService> fs;
   mtree = fs->make<TTree>("Events", "Events");
  
  mtree->Branch("numberelectron",&numelectron);   
  mtree->GetBranch("numberelectron")->SetTitle("number of electrons");
  mtree->Branch("electron_e",&electron_e);
  mtree->GetBranch("electron_e")->SetTitle("electron energy");
  mtree->Branch("electron_pt",&electron_pt);
  mtree->GetBranch("electron_pt")->SetTitle("electron transverse momentum");
  mtree->Branch("electron_px",&electron_px);
  mtree->GetBranch("electron_px")->SetTitle("electron momentum x-component");
  mtree->Branch("electron_py",&electron_py);
  mtree->GetBranch("electron_py")->SetTitle("electron momentum y-component");
  mtree->Branch("electron_pz",&electron_pz);
  mtree->GetBranch("electron_pz")->SetTitle("electron momentum z-component");
  mtree->Branch("electron_eta",&electron_eta);
  mtree->GetBranch("electron_eta")->SetTitle("electron pseudorapidity");
  mtree->Branch("electron_phi",&electron_phi);
  mtree->GetBranch("electron_phi")->SetTitle("electron polar angle");
~~~
{: .language-cpp}

Essentially all of the EDAnalyzers clear the variable containers and loop over the objects present in the event.  

Finally, many of the most important kinematic quantities defining a physics object are accessed in a common way across all the objects. 
Most objects have associated energy-momentum vectors, typically constructed using **transverse momentum, pseudorapdity, azimuthal angle, and mass or energy**.

~~~
   numelectron = 0;
   electron_e.clear();
   electron_pt.clear();
   electron_px.clear();
   electron_py.clear();
   electron_pz.clear();
   electron_eta.clear();
   electron_phi.clear();
   electron_ch.clear();
   electron_iso.clear();
   electron_veto.clear();//
   electron_isLoose.clear();
   electron_isMedium.clear();
   electron_isTight.clear();
   electron_dxy.clear();
   electron_dz.clear();
   electron_dxyError.clear();
   electron_dzError.clear();
   electron_ismvaLoose.clear();
   electron_ismvaTight.clear();

    for (const pat::Electron &el : *electrons)
    {
      electron_e.push_back(el.energy());
      electron_pt.push_back(el.pt());
      electron_px.push_back(el.px());
      electron_py.push_back(el.py());
      electron_pz.push_back(el.pz());
      electron_eta.push_back(el.eta());
      electron_phi.push_back(el.phi());
~~~
{: .language-cpp}

{% include links.md %}

./rootcom postAnalyzer_mc_JESPES postAnalyzer_mc
./rootcom postAnalyzer_mc_Wnugamma_JESPES postAnalyzer_Wnugamma
#./rootcom postAnalyzer_mc_ZLLGamma_JESPES postAnalyzer_ZLLGamma
#./rootcom postAnalyzer_mc_ZnunuGamma_JESPES_LO postAnalyzer_ZnunuGamma
#========Below is the path at KISTI 

./MakeCondorFiles.sh postAnalyzer_mc root://cms-t2-se01.sdfarm.kr:1096//store/user/sdogra/Summer20ULMiniAOD_Dec2021/2017_v2/Summer20UL2017MiniAODv2/TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/crab_job_TTGJets/220317_010914/0000/  TTGJets_new.root -1 1000 log_TTGJets ewk_corr.root egammaEffi.txt_EGM2D_uncer_mySF.root

./MakeCondorFiles.sh postAnalyzer_mc root://cms-t2-se01.sdfarm.kr:1096//store/user/sdogra/Summer20ULMiniAOD_Dec2021/2017_v2/Summer20UL2017MiniAODv2/TGJets_TuneCP5_13TeV-amcatnlo-madspin-pythia8/crab_job_TGJets/220317_010930/0000/ TGJets_new.root -1 1000 log_TGJets ewk_corr.root egammaEffi.txt_EGM2D_uncer_mySF.root

./MakeCondorFiles.sh postAnalyzer_mc root://cms-t2-se01.sdfarm.kr:1096//store/user/sdogra/Summer20ULMiniAOD_Dec2021/2017_v2/Summer20UL2017MiniAODv2/DiPhotonJets_MGG-80toInf_TuneCP5_13TeV-amcatnloFXFX-pythia8/crab_job_DiPhotonJets/220317_011208/0000/ DiPhotonJets_new.root -1 1000 log_DiPhotonJets ewk_corr.root egammaEffi.txt_EGM2D_uncer_mySF.root

./MakeCondorFiles.sh postAnalyzer_mc root://cms-t2-se01.sdfarm.kr:1096//store/user/sdogra/Summer20ULMiniAOD_Dec2021/2017_v2/Summer20UL2017MiniAODv2/WW_TuneCP5_13TeV-pythia8/crab_job_WW/220317_010944/0000/ WW_new.root -1 1000 log_WW ewk_corr.root egammaEffi.txt_EGM2D_uncer_mySF.root

./MakeCondorFiles.sh postAnalyzer_mc root://cms-t2-se01.sdfarm.kr:1096//store/user/sdogra/Summer20ULMiniAOD_Dec2021/2017_v2/Summer20UL2017MiniAODv2/WZ_TuneCP5_13TeV-pythia8/crab_job_WZ/220317_010959/0000/ WZ_new.root -1 1000 log_WZ ewk_corr.root egammaEffi.txt_EGM2D_uncer_mySF.root

./MakeCondorFiles.sh postAnalyzer_mc root://cms-t2-se01.sdfarm.kr:1096//store/user/sdogra/Summer20ULMiniAOD_Dec2021/2017_v2/Summer20UL2017MiniAODv2/ZZ_TuneCP5_13TeV-pythia8/crab_job_ZZ/220317_011013/0000/ ZZ_new.root -1 1000 log_ZZ ewk_corr.root egammaEffi.txt_EGM2D_uncer_mySF.root

./MakeCondorFiles.sh postAnalyzer_mc root://cms-t2-se01.sdfarm.kr:1096//store/user/sdogra/Summer20ULMiniAOD_Dec2021/2017_v2/Summer20UL2017MiniAODv2/WToMuNu_M-100_TuneCP5_13TeV-pythia8/crab_job_WToMuNu_M-100/220317_010845/0000/ WToMuNu_M.root -1 1000 log_WToMuNu_M ewk_corr.root egammaEffi.txt_EGM2D_uncer_mySF.root

./MakeCondorFiles.sh postAnalyzer_mc root://cms-t2-se01.sdfarm.kr:1096//store/user/sdogra/Summer20ULMiniAOD_Dec2021/2017_v2/Summer20UL2017MiniAODv2/WToTauNu_M-100_TuneCP5_13TeV-pythia8-tauola/crab_job_WToTauNu_M-100/220317_010900/0000/ WToTauNu_new.root -1 1000 log_WToTauNu ewk_corr.root egammaEffi.txt_EGM2D_uncer_mySF.root


./MakeCondorFiles.sh postAnalyzer_mc root://cms-t2-se01.sdfarm.kr:1096//store/user/sdogra/Summer20ULMiniAOD_Dec2021/2017_v2/Summer20UL2017MiniAODv2/GJets_DR-0p4_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/crab_job_GJets_DR-0p4_HT-100To200/220312_185640/0000/ GJets_DR-0p4_HT-100To200_new.root -1 10000 log_GJets_DR-0p4_HT-100To200 ewk_corr.root egammaEffi.txt_EGM2D_uncer_mySF.root

./MakeCondorFiles.sh postAnalyzer_mc root://cms-t2-se01.sdfarm.kr:1096//store/user/sdogra/Summer20ULMiniAOD_Dec2021/2017_v2/Summer20UL2017MiniAODv2/GJets_DR-0p4_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/crab_job_GJets_DR-0p4_HT-200To400/220310_124510/0000/ GJets_DR-0p4_HT-200To400_new.root -1 10000  log_GJets_DR-0p4_HT-200To400 ewk_corr.root egammaEffi.txt_EGM2D_uncer_mySF.root

./MakeCondorFiles.sh postAnalyzer_mc root://cms-t2-se01.sdfarm.kr:1096//store/user/sdogra/Summer20ULMiniAOD_Dec2021/2017_v2/Summer20UL2017MiniAODv2/GJets_DR-0p4_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/crab_job_GJets_DR-0p4_HT-400To600/220224_072114/0000/ GJets_DR-0p4_HT-400To600_new.root -1 10000  log_GJets_DR-0p4_HT-400To600 ewk_corr.root egammaEffi.txt_EGM2D_uncer_mySF.root

./MakeCondorFiles.sh postAnalyzer_mc root://cms-t2-se01.sdfarm.kr:1096//store/user/sdogra/Summer20ULMiniAOD_Dec2021/2017_v2/Summer20UL2017MiniAODv2/GJets_DR-0p4_HT-600ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/crab_job_GJets_DR-0p4_HT-600ToInf/220224_072136/0000/ GJets_DR-0p4_HT-600ToInf_new.root  -1 10000  log_GJets_DR-0p4_HT-600ToInf ewk_corr.root egammaEffi.txt_EGM2D_uncer_mySF.root

./MakeCondorFiles.sh postAnalyzer_mc root://cms-t2-se01.sdfarm.kr:1096//store/user/sdogra/Summer20ULMiniAOD_Dec2021/2017_v2/Summer20UL2017MiniAODv2/GJets_HT-40To100_TuneCP5_13TeV-madgraphMLM-pythia8/crab_job_GJets_DR-0p4_HT-40To100/220312_185624/0000/ GJets_DR-0p4_HT-40To100_new.root -1 10000 log_GJets_DR-0p4_HT-40To100 ewk_corr.root egammaEffi.txt_EGM2D_uncer_mySF.root

#./MakeCondorFiles.sh postAnalyzer_ZLLGamma root://cms-t2-se01.sdfarm.kr:1096//store/user/sdogra/Summer20ULMiniAOD_Dec2021/2017_v2/Summer20UL2017MiniAODv2/ZLLGJets_MonoPhoton_PtG-130_TuneCP5_13TeV-amcatnloFXFX-pythia8/crab_job_ZLLGJets/220317_010802/0000/ ZLLGJets.root -1 1000 log_ZLLGJets  ewk_corr.root egammaEffi.txt_EGM2D_uncer_mySF.root

./MakeCondorFiles.sh postAnalyzer_Wnugamma root://cms-t2-se01.sdfarm.kr:1096//store/user/sdogra/Summer20ULMiniAOD_Dec2021/2017_v2/Summer20UL2017MiniAODv2/WGJets_MonoPhoton_PtG-130_TuneCP5_13TeV-madgraph-pythia8/crab_job_WGJets_MonoPhoton_PtG_130/220317_010816/0000/ WGJets_new.root -1 1000 log_WGJets  ewk_corr.root egammaEffi.txt_EGM2D_uncer_mySF.root


#./postAnalyzer_ZLLGamma /T2_KR_KISTI/store/user/sdogra/Summer20ULMiniAOD_Dec2021/2017_v2/Summer20UL2017MiniAODv2/ZLLGJets_MonoPhoton_PtG-130_TuneCP5_13TeV-amcatnloFXFX-pythia8/crab_job_ZLLGJets/220317_010802/0000/ ZLLGJets.root -1 1000 > log_ZLLGJets  


#./postAnalyzer_Wnugamma /T2_KR_KISTI/store/user/sdogra/Summer20ULMiniAOD_Dec2021/2017_v2/Summer20UL2017MiniAODv2/WGJets_MonoPhoton_PtG-130_TuneCP5_13TeV-madgraph-pythia8/crab_job_WGJets_MonoPhoton_PtG_130/220317_010816/0000/ WGJets.root -1 1000 > log_WGJets



#./postAnalyzer_ZnunuGamma /cms/scratch/sdogra/mPhoton/Znunu_LO/Znunu17/ Znunugamma17_LO.root -1 1000  > log_Znunugamma17




# Helper functions for differential expression analyses

# Each function iterates through the tiered cell type classification and in each cell type at each tier
 # finds the genes differentially expressed between EE and the specified control cohorts

# Finds genes differentially expressed between EE and all three control cohorts
find_markers_three = function(sobj){
  Idents(sobj) <- "study" # I think that for now doing DE by study is probably still the best thing to be doing, based on healthy being similar to active all
  tones = unique(sobj$tier1)
  print(tones)
  for(tone in tones){
    print(paste0("tier1: ",tone))
    print(typeof(tone))
    stone = sobj[,colnames(sobj)[sobj$tier1==tone]]
    if("EE" %in% stone$study & min(table(stone$study)) > 2 & length(unique(stone$study))>1){
      markers_all = FindMarkers(stone,logfc.threshold = 0.1,min.pct = 0.025, ident.1="EE")
      write.csv(markers_all,paste0("tier1",tone,"all.csv"))
    }
    ttwos = unique(sobj$tier2[sobj$tier1==tone]) # ok now just looking at the subtree
    if(length(ttwos)>1){
      for(ttwo in ttwos){
        print(paste0("tier2: ",ttwo))
        sttwo = stone[,colnames(stone)[stone$tier2==ttwo]]
        if("EE" %in% sttwo$study & min(table(sttwo$study)) > 2 & length(unique(sttwo$study))>1){
          markers_all = FindMarkers(sttwo,logfc.threshold = 0.1,min.pct = 0.025, ident.1="EE")
          write.csv(markers_all,paste0("tier1",tone,"_tier2",ttwo,"all.csv"))
        }
        tthrees = unique(sobj$tier3[sobj$tier1==tone & sobj$tier2==ttwo])
        if(length(tthrees)>1){
          for(tthree in tthrees){
            print(paste0("tier3: ",tthree))
            stthree = sttwo[,colnames(sttwo)[sttwo$tier3==tthree]]
            if("EE" %in% stthree$study & min(table(stthree$study)) > 2  & length(unique(stthree$study))>1){
              markers_all = FindMarkers(stthree,logfc.threshold = 0.1,min.pct = 0.025, ident.1="EE")
              write.csv(markers_all,paste0("tier1",tone,"_tier2",ttwo,"_tier3",tthree,"all.csv"))
            }
            tfours = unique(sobj$tier4[sobj$tier1==tone & sobj$tier2==ttwo & sobj$tier3==tthree])
            if(length(tfours)>1){
              for(tfour in tfours){
                print(paste0("tier4: ",tfour))
                stfour = stthree[,colnames(stthree)[stthree$tier4==tfour]]
                if("EE" %in% stfour$study & min(table(stfour$study)) > 2 & length(unique(stfour$study))>1){
                  markers_all = FindMarkers(stfour,logfc.threshold = 0.1,min.pct = 0.025, ident.1="EE")
                  write.csv(markers_all,paste0("tier1",tone,"_tier2",ttwo,"_tier3",tthree,"_tier4",tfour,"all.csv"))
                }
                
              }
            }
          }
          
        }
      }
    }
  }
}

# Finds genes differentially expressed between EE and both U.S. control cohorts
find_markers_two = function(sobj){
  Idents(sobj) <- "study"  
  tones = unique(sobj$tier1)
  print(tones)
  for(tone in tones){
    print(paste0("tier1: ",tone))
    print(typeof(tone))
    stone = sobj[,colnames(sobj)[sobj$tier1==tone]]
    if("EE" %in% stone$study & min(table(stone$study)) > 2 & length(unique(stone$study))>1){
      markers_both = FindMarkers(stone,logfc.threshold = 0.1,min.pct = 0.025, ident.1="EE")
      write.csv(markers_both,paste0("tier1",tone,"both.csv"))
    }
    ttwos = unique(sobj$tier2[sobj$tier1==tone]) 
    if(length(ttwos)>1){
      for(ttwo in ttwos){
        print(paste0("tier2: ",ttwo))
        sttwo = stone[,colnames(stone)[stone$tier2==ttwo]]
        if("EE" %in% sttwo$study & min(table(sttwo$study)) > 2 & length(unique(sttwo$study))>1){
          markers_both = FindMarkers(sttwo,logfc.threshold = 0.1,min.pct = 0.025, ident.1="EE")
          write.csv(markers_both,paste0("tier1",tone,"_tier2",ttwo,"both.csv"))
        }
        tthrees = unique(sobj$tier3[sobj$tier1==tone & sobj$tier2==ttwo])
        if(length(tthrees)>1){
          for(tthree in tthrees){
            print(paste0("tier3: ",tthree))
            stthree = sttwo[,colnames(sttwo)[sttwo$tier3==tthree]]
            if("EE" %in% stthree$study & min(table(stthree$study)) > 2  & length(unique(stthree$study))>1){
              markers_both = FindMarkers(stthree,logfc.threshold = 0.1,min.pct = 0.025, ident.1="EE")
              write.csv(markers_both,paste0("tier1",tone,"_tier2",ttwo,"_tier3",tthree,"both.csv"))
            }
            tfours = unique(sobj$tier4[sobj$tier1==tone & sobj$tier2==ttwo & sobj$tier3==tthree])
            if(length(tfours)>1){
              for(tfour in tfours){
                print(paste0("tier4: ",tfour))
                stfour = stthree[,colnames(stthree)[stthree$tier4==tfour]]
                if("EE" %in% stfour$study & min(table(stfour$study)) > 2 & length(unique(stfour$study))>1){
                  markers_both = FindMarkers(stfour,logfc.threshold = 0.1,min.pct = 0.025, ident.1="EE")
                  write.csv(markers_both,paste0("tier1",tone,"_tier2",ttwo,"_tier3",tthree,"_tier4",tfour,"both.csv"))
                }
                
              }
            }
          }
          
        }
      }
    }
  }
}

# Finds genes differentially expressed between EE and Durban cohort
find_markers_Durban = function(sobj){
  Idents(sobj) <- "study" 
  tones = unique(sobj$tier1)
  print(tones)
  for(tone in tones){
    print(paste0("tier1: ",tone))
    print(typeof(tone))
    
    stone = sobj[,colnames(sobj)[sobj$tier1==tone]]
    
    if(length(unique(stone$study))>1 & min(table(stone$study)) > 2){
      markers_Durban = FindMarkers(stone,logfc.threshold = 0.1,min.pct = 0.025, ident.1="EE",ident.2="Durban")
      write.csv(markers_Durban,paste0("tier1",tone,"Durban.csv"))
    }
    ttwos = unique(sobj$tier2[sobj$tier1==tone])
    if(length(ttwos)>1){
      for(ttwo in ttwos){
        print(paste0("tier2: ",ttwo))
        sttwo = stone[,colnames(stone)[stone$tier2==ttwo]]
        if(length(unique(sttwo$study))>1 & min(table(sttwo$study)) > 2){
          markers_Durban = FindMarkers(sttwo,logfc.threshold = 0.1,min.pct = 0.025, ident.1="EE",ident.2="Durban")
          write.csv(markers_Durban,paste0("tier1",tone,"_tier2",ttwo,"Durban.csv"))
        }
        tthrees = unique(sobj$tier3[sobj$tier1==tone & sobj$tier2==ttwo])
        if(length(tthrees)>1){
          for(tthree in tthrees){
            print(paste0("tier3: ",tthree))
            stthree = sttwo[,colnames(sttwo)[sttwo$tier3==tthree]]
            if(length(unique(stthree$study))>1 & min(table(stthree$study)) > 2){
              markers_Durban = FindMarkers(stthree,logfc.threshold = 0.1,min.pct = 0.025, ident.1="EE",ident.2="Durban")
              write.csv(markers_Durban,paste0("tier1",tone,"_tier2",ttwo,"_tier3",tthree,"Durban.csv"))
            }
            tfours = unique(sobj$tier4[sobj$tier1==tone & sobj$tier2==ttwo & sobj$tier3==tthree])
            if(length(tfours)>1){
              for(tfour in tfours){
                print(paste0("tier4: ",tfour))
                stfour = stthree[,colnames(stthree)[stthree$tier4==tfour]]
                if(length(unique(stfour$study))>1 & min(table(stfour$study)) > 2){
                  markers_Durban = FindMarkers(stfour,logfc.threshold = 0.1,min.pct = 0.025, ident.1="EE",ident.2="Durban")
                  write.csv(markers_Durban,paste0("tier1",tone,"_tier2",ttwo,"_tier3",tthree,"_tier4",tfour,"Durban.csv"))
                }
              }
            }
          }
          
        }
      }
    }
  }
}

# Finds genes differentially expressed between EE and EoE cohort
find_markers_EoE = function(sobj){
  Idents(sobj) <- "study" 
  tones = unique(sobj$tier1)
  print(tones)
  for(tone in tones){
    print(paste0("tier1: ",tone))
    print(typeof(tone))
    
    stone = sobj[,colnames(sobj)[sobj$tier1==tone]]
    
    if(length(unique(stone$study))>1 & min(table(stone$study)) > 2){
      markers_EoE = FindMarkers(stone,logfc.threshold = 0.1,min.pct = 0.025, ident.1="EE",ident.2="EoE")
      write.csv(markers_EoE,paste0("tier1",tone,"EoE.csv"))
    }
    ttwos = unique(sobj$tier2[sobj$tier1==tone]) 
    if(length(ttwos)>1){
      for(ttwo in ttwos){
        print(paste0("tier2: ",ttwo))
        sttwo = stone[,colnames(stone)[stone$tier2==ttwo]]
        if(length(unique(sttwo$study))>1 & min(table(sttwo$study)) > 2){
          markers_EoE = FindMarkers(sttwo,logfc.threshold = 0.1,min.pct = 0.025, ident.1="EE",ident.2="EoE")
          write.csv(markers_EoE,paste0("tier1",tone,"_tier2",ttwo,"EoE.csv"))
        }
        tthrees = unique(sobj$tier3[sobj$tier1==tone & sobj$tier2==ttwo])
        if(length(tthrees)>1){
          for(tthree in tthrees){
            print(paste0("tier3: ",tthree))
            stthree = sttwo[,colnames(sttwo)[sttwo$tier3==tthree]]
            if(length(unique(stthree$study))>1 & min(table(stthree$study)) > 2){
              markers_EoE = FindMarkers(stthree,logfc.threshold = 0.1,min.pct = 0.025, ident.1="EE",ident.2="EoE")
              write.csv(markers_EoE,paste0("tier1",tone,"_tier2",ttwo,"_tier3",tthree,"EoE.csv"))
            }
            tfours = unique(sobj$tier4[sobj$tier1==tone & sobj$tier2==ttwo & sobj$tier3==tthree])
            if(length(tfours)>1){
              for(tfour in tfours){
                print(paste0("tier4: ",tfour))
                stfour = stthree[,colnames(stthree)[stthree$tier4==tfour]]
                if(length(unique(stfour$study))>1 & min(table(stfour$study)) > 2){
                  markers_EoE = FindMarkers(stfour,logfc.threshold = 0.1,min.pct = 0.025, ident.1="EE",ident.2="EoE")
                  write.csv(markers_EoE,paste0("tier1",tone,"_tier2",ttwo,"_tier3",tthree,"_tier4",tfour,"EoE.csv"))
                }
              }
            }
          }
          
        }
      }
    }
  }
}

# Finds genes differentially expressed between EE and Resection cohort
find_markers_Resection = function(sobj){
  Idents(sobj) <- "study"
  tones = unique(sobj$tier1)
  print(tones)
  for(tone in tones){
    print(paste0("tier1: ",tone))
    print(typeof(tone))
    
    stone = sobj[,colnames(sobj)[sobj$tier1==tone]]
    
    if(length(unique(stone$study))>1 & min(table(stone$study)) > 2){
      markers_Resection = FindMarkers(stone,logfc.threshold = 0.1,min.pct = 0.025, ident.1="EE",ident.2="Resection")
      write.csv(markers_Resection,paste0("tier1",tone,"Resection.csv"))
    }
    ttwos = unique(sobj$tier2[sobj$tier1==tone]) 
    if(length(ttwos)>1){
      for(ttwo in ttwos){
        print(paste0("tier2: ",ttwo))
        sttwo = stone[,colnames(stone)[stone$tier2==ttwo]]
        if(length(unique(sttwo$study))>1 & min(table(sttwo$study)) > 2){
          markers_Resection = FindMarkers(sttwo,logfc.threshold = 0.1,min.pct = 0.025, ident.1="EE",ident.2="Resection")
          write.csv(markers_Resection,paste0("tier1",tone,"_tier2",ttwo,"Resection.csv"))
        }
        tthrees = unique(sobj$tier3[sobj$tier1==tone & sobj$tier2==ttwo])
        if(length(tthrees)>1){
          for(tthree in tthrees){
            print(paste0("tier3: ",tthree))
            stthree = sttwo[,colnames(sttwo)[sttwo$tier3==tthree]]
            if(length(unique(stthree$study))>1 & min(table(stthree$study)) > 2){
              markers_Resection = FindMarkers(stthree,logfc.threshold = 0.1,min.pct = 0.025, ident.1="EE",ident.2="Resection")
              write.csv(markers_Resection,paste0("tier1",tone,"_tier2",ttwo,"_tier3",tthree,"Resection.csv"))
            }
            tfours = unique(sobj$tier4[sobj$tier1==tone & sobj$tier2==ttwo & sobj$tier3==tthree])
            if(length(tfours)>1){
              for(tfour in tfours){
                print(paste0("tier4: ",tfour))
                stfour = stthree[,colnames(stthree)[stthree$tier4==tfour]]
                if(length(unique(stfour$study))>1 & min(table(stfour$study)) > 2){
                  markers_Resection = FindMarkers(stfour,logfc.threshold = 0.1,min.pct = 0.025, ident.1="EE",ident.2="Resection")
                  write.csv(markers_Resection,paste0("tier1",tone,"_tier2",ttwo,"_tier3",tthree,"_tier4",tfour,"Resection.csv"))
                }
              }
            }
          }
          
        }
      }
    }
  }
}
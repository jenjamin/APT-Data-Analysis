# This function takes the IVAS cluster search csv output
# It separates ions into their constituent elements
# For RPV steel users, all Fm ions within the cluster are assumed to be Ni
# It then converts the cluster information into a new dataframe called "ClusterAtomCount"
# Written by Ben Jenkins Dec 2018

CSVInputConverter <- function(DataframeName){
  ClusterInformation <- DataframeName %>%
    select("X",contains("Ions"), contains("Count")) %>%
    rename_(.dots=setNames(names(.), (gsub(".Count", "", names(.)))))
  
  #### Creates columns for elements if they are not present. Column values will equal 0 ####
  fncols <- function(data, cname) {
    add <-cname[!cname%in%names(data)]
    
    if(length(add)!=0) data[add] <- 0
    ClusterInfoAllElements <<- data
  }
  fncols(ClusterInformation, c("Fm", "Cu", "Ni", "Mn", "Si", "P"))
  
  #### Converts ion counts to individual atom counts ####
  require(PeriodicTable)
  require(InterpretMSSpectrum)
  data("periodicTable")
  ElementsList <- c(unique(periodicTable$symb))
  ClusterInformationConverted <- ClusterInfoAllElements
  for(Ions in colnames(ClusterInfoAllElements %>% select(-"X",-contains("Ions")))){
    if(Ions %in% unique(periodicTable$symb) == FALSE){
      SplitIon <- CountChemicalElements(Ions)
      i <- 1
      for(DecomposedElements in names(SplitIon[])){
        ClusterInformationConverted <- ClusterInformationConverted %>%
          mutate(!!paste((DecomposedElements), sep = "_") := (!!rlang::sym(paste((DecomposedElements), sep = "_"))) + 
                   ((!!rlang::sym(paste((Ions), sep = "_"))) * as.numeric(SplitIon[i]))) 
        i < i + 1 
      }
      ClusterInformationConverted <<- ClusterInformationConverted %>% select(-Ions)
    }
  }
  #rm(DecomposedElements, i, Ions, SplitIon, ElementsList) # Remove vars from global envir
  

  ClusterInformationConvertedCount <<- ClusterInformationConverted %>%
    mutate(TotalAtoms = rowSums(select(.,-X, -Solute.Ions, -Ranged.Ions, -Total.Ions)),
           Ni = Ni + Fm) %>% # Assumes all Fm in clusters in Ni
    select(-Fm) %>%
    mutate(Other = TotalAtoms - (Fe + Si + Cu + Ni + Mn + P),
           NonFeAtoms = TotalAtoms-Fe) %>%
    select(X, Solute.Ions, Ranged.Ions, TotalAtoms, NonFeAtoms,
           Fe, Si, Cu, Ni, Mn, P, Other)
  
  #### Create dataframe ready to plot bar charts ####  
  ClusterInfoForBarChart <<- ClusterInformationConverted %>%
    mutate(TotalAtoms = rowSums(select(.,-X, -Solute.Ions, -Ranged.Ions, -Total.Ions)),
           Ni = Ni + Fm) %>% # Assumes all Fm in clusters in Ni
    select(-Fm) %>%
    mutate(Other = TotalAtoms - (Fe + Si + Cu + Ni + Mn + P),
           NonFeAtoms = TotalAtoms-Fe) %>%
    select(X, Solute.Ions, Ranged.Ions, TotalAtoms, NonFeAtoms,
           Fe, Si, Cu, Ni, Mn, P, Other)%>%
    filter(grepl("Cluster",X)) %>%
    mutate(Fe = 100 * Fe/Ranged.Ions,
           Si = 100 * Si/NonFeAtoms, 
           Cu = 100 * Cu/NonFeAtoms, 
           Ni = 100 * Ni/NonFeAtoms, 
           Mn = 100 * Mn/NonFeAtoms, 
           P = 100 * P/NonFeAtoms, 
           Other = 100 * Other/NonFeAtoms,
           TotalAtoms = TotalAtoms + (parse_number(as.character(X))/10000)) %>%
    gather(Element, `Concentration (at.%)`, Fe, Si, Cu, Ni, Mn, P, Other)
  
}
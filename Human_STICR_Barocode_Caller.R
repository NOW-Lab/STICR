library(data.table)
library(parallel)
library(Matrix)
library(ggplot2)
library(rlist)
######Functions######
VBC_Jaccard_Precompiled_list_function<-function(x){
  index.set<-VBC_Barcode_Index.comb[x,]
  Intersection.res<-length(intersect(Cell.list[[index.set$X1]],Cell.list[[index.set$X2]]))
  Union.res<-length(union(Cell.list[[index.set$X1]],Cell.list[[index.set$X2]]))
  Jaccard.sim<-Intersection.res/Union.res
  Jacc.Union.res<-paste(Jaccard.sim,Union.res,sep="_")
  return(Jacc.Union.res)
}

Confirmation_function<-function(x){
  working.set<-as.numeric(High_Confidence_Superinfection[x,1:2])
  Quote.part1<-dQuote(paste("\\b",working.set[1],"\\b", sep=""))
  Quote.part2<-substr(Quote.part1,2,nchar(Quote.part1)-1)
  returns<-grep(Quote.part2,SuperInfection.ls)
  confirmation<-ifelse(working.set[2] %in% unlist(SuperInfection.ls[returns]), 1,0)
  return(confirmation)
} #confirms that each pair is present in final superinfection set

Super_Infection_Matching_Function<-function(x){
  working.set<-SuperInfection.ls[[x]]
  VBCs<-as.character(VBC_Barcode_Index$VBC[VBC_Barcode_Index$Index %in% working.set])
  return(VBCs)
}

Barcode_Collapse_Function<-function(x){
  VBCs<-as.character(working.cell.split[[x]]$barcode)
  return(VBCs)
}

Barcode_Calling_Function<-function(x){
  working.set<-working.cell.split[[x]]
  barcode.set<-as.character(working.set$barcode)
  hits<-list()
  for (i in barcode.set){
    hits<-list.append(hits,grep(i,Full_VBC.Single.ls, fixed = TRUE))
  }
  result<-ifelse(length(unique(unlist(hits)))>1,0, unique(unlist(hits)) )
  return(result)
} #Multiplet result will return a value of -1, Single Cell Clones will return a value of -2

Multi_VBC_count_function<-function(x){
  working.set<-working.cell.split[[x]]
  barcode.set<-as.character(working.set$barcode)
  VBC_count<-length(unlist(barcode.set))
  return(VBC_count)
} #checking to see what % of cells belonging to multi-VBC clones contain all VBCs

######Import and Format Data for Analysis######

#Import VBC calls
Full_BC.df<-readRDS("/STICR_Barcode_Extractor_outputFolder/Final_Barcodes.tsv") 

#Import "Final_Barcodes.tsv" file generated from the STICR Barcode Extractor Script as a data.frame named "Full_BC.df"
#Add an additional column called "Library" that describes the 10X transcriptomic library this STICR Barcode data was generated from. Value should be a character string.
#Dataframes from multiple libraries can be concatenated together into data.frame named "Full_BC.df".
#Final data.frame should contain 4 columns with names: "CBC","barcode","UMI_Count","Library"


#Create data.frame of "valid" CBCs- i.e. those in scRNA-seq transcriptomic set passing QC cutoffs
Clone_Seq_BC.df<-#valid CBC set
#Clone_Seq_BC.df should have 2 columns titled: "Library" and "CBC"
#Library values should be character strings and match values in Full_BC.df$Library
#CBC values should be character strings and match formating used in Full_BC.df$CBC. Of note, depending on how datasets were processed/merged in Seurat, CBC values might have a "-N" appended to the end that you will want to match.

#Subset VBC dataframe by CBCs detected in Transcriptomic set 
Full_BC.df<-Full_BC.df[Full_BC.df$CBC %in% Clone_Seq_BC.df$CBC,]

#Retain CBC/STICR combinations with at least 5 UMI count
Full_BC.df<-Full_BC.df[Full_BC.df$UMI_Count>4,] 



#Set Variables
numCores<-10
library.set<-unique(Full_BC.df$Library_Full) #Create vector of library names

#Analysis

#Checking for Superinfection
for (i in library.set){
  working.lib<-data.frame(Full_BC.df[Full_BC.df$Library_Full == i,]) #Subset VBC Count Matrix by Library
  VBC_Barcode_Index<-data.frame(VBC=sort(unique(working.lib[,2])),Index=1:length(unique(working.lib[,2]))) #Create a numerical VBC index
  VBC_Single_Cell.vec<-names(table(working.lib$barcode)[table(working.lib$barcode)==1]) #create character vector of VBCs that appear in a single cell
  VBC_Barcode_Index$Single_Cell_Clone<-ifelse(VBC_Barcode_Index$VBC %in% VBC_Single_Cell.vec,1,0) #Add a metadata column to numerical VBC index indicating whether a VBC appears in just one cell (value=1) or not (value=0)
  VBC_Barcode_Index.comb<-data.frame(matrix(unlist(combn(VBC_Barcode_Index$Index[VBC_Barcode_Index$Single_Cell_Clone==0],2,simplify=FALSE)),ncol=2, byrow=TRUE)) #make all unique pairwise between non-singlet VBCs
  Cell.list<-split(working.lib[,1],working.lib[,2]) #split data.frame by VBC
  VBC_Barcode_Index.comb$jaccard<-mclapply(1:dim(VBC_Barcode_Index.comb)[1], VBC_Jaccard_Precompiled_list_function, mc.cores = numCores) #Compute jaccard smilarity index for each VBC pair (defined here as the Union/Intersection of cells with VBCs being tested)
  VBC_Barcode_Index.comb$Union<- as.numeric(sub('.*_', '', VBC_Barcode_Index.comb$jaccard)) #extract and convert the size of Union for the VBC pair from the contatenated character output from previous function to a numeric vector
  VBC_Barcode_Index.comb$jaccard<-as.numeric(sub('_.*', '', VBC_Barcode_Index.comb$jaccard)) #extract and convert the size of Jaccard similarity index for the VBC pair from the contatenated character output from previous function to a numeric vector
  VBC_Jaccard.single<-ggplot(VBC_Barcode_Index.comb[VBC_Barcode_Index.comb$jaccard > 0 & VBC_Barcode_Index.comb$Union > 9,], aes(jaccard)) + geom_histogram(binwidth = 0.05, color="black", fill="#4575b4") + scale_y_continuous(trans = "log2") + geom_vline(xintercept = 0.55, linetype="dotted", color = "red", size=1.5) + xlim(0,1) #make a ggplot2 histogram of jaccard similarity index frequency, subsetting for VBC pairs whose Union value is >4 and Jaccard Similarity Index is >0.
  title<-paste0("Jaccard_Histogram_",i,".pdf") #Create a title for the histogram complete with file extension (.pdf)
  pdf(title, width = 4, height = 4) 
  print(VBC_Jaccard.single + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + ggtitle(title) + xlab("Jaccard Similarity") + xlab("Log2(Count)"))
  dev.off()
}

#Check the histogram to determine if you want to change the cutoff for including a superinfection. 0.55 was chosen for our analysis as our datasets as this threshold appeared to indicate the start a second "peak" in the histogram data.
#From here, users can loop through the rest of this code for each library in the Full_BC.df data.frame


High_Confidence_Superinfection<-VBC_Barcode_Index.comb[VBC_Barcode_Index.comb$jaccard >= 0.55 &VBC_Barcode_Index.comb$Union>9,1:2] #Create set of VBC pairs that have a Jaccard Similarity Index higher than a pre-specified cutoff and that occur with a Union of 10 or more cells.
starting.vector<-c(High_Confidence_Superinfection$X1,High_Confidence_Superinfection$X2) #create a list of all VBCs suspected to be in superinfections
SuperInfection.ls<-list()
while(length(starting.vector)>0) {
  t=1
  working.set<-starting.vector[1]
  while (t!=0) {
    start.set<-working.set
    X1.set<-High_Confidence_Superinfection$X1[High_Confidence_Superinfection$X2 %in% start.set]
    X2.set<-High_Confidence_Superinfection$X2[High_Confidence_Superinfection$X1 %in% start.set]
    end.set<-unique(c(start.set,X1.set,X2.set))
    starting.vector<-setdiff(starting.vector,end.set)
    t<-length(end.set)-length(start.set)
    if (t==0){
      SuperInfection.ls<-list.append(SuperInfection.ls,end.set)
    }
    working.set<-end.set
  }
} #Determines whether suspected VBC superinfection pairs have additional overlap (ie. 3+ VBCs per superinfection)

if(length(SuperInfection.ls)>0){
  High_Confidence_Superinfection$confirmation<-as.numeric(lapply(1:dim(High_Confidence_Superinfection)[1],Confirmation_function )) #table of the High_Confidence_Superinfection$confirmation value should produce all 1's
  Superinfection_Grouping.score.df<-data.frame(table(High_Confidence_Superinfection$confirmation))
  Superinfection_Grouping.score<-(Superinfection_Grouping.score.df$Freq[Superinfection_Grouping.score.df$Var1==1])
  SuperInfection.VBCs.ls<-lapply(1:length(SuperInfection.ls),Super_Infection_Matching_Function) #Replaces Index with actual barcodes
  Superinfection_SingleCellClone_check<-intersect(VBC_Single_Cell.vec,unlist(SuperInfection.VBCs.ls)) #makes sure that VBCs found in single cell not included in Superinfection
}



VBC_Barcode_Index$SuperInfection<-ifelse(VBC_Barcode_Index$Index %in% unlist(SuperInfection.ls),1,0) #Add a metadata column to numerical VBC index indicating whether a VBC appears as part of a superinfection (value=1) or not (value=0)
MOI_1_VBC<-as.character(VBC_Barcode_Index$VBC[VBC_Barcode_Index$SuperInfection==0 & VBC_Barcode_Index$Single_Cell_Clone==0]) #Imputed Single VBC clones with >1 cell
working.cell.split<-split(working.lib,working.lib$CBC) #split data.frame by CBC
CBC_VBC_Calls<-lapply(names(working.cell.split),Barcode_Collapse_Function)
names(CBC_VBC_Calls)<-names(working.cell.split)
Full_VBC.ls<-list.append(MOI_1_VBC,SuperInfection.VBCs.ls) #full list of valid VBCs
Full_VBC.Single.ls<-c(Full_VBC.ls,VBC_Single_Cell.vec)
Barcodes.aligned<-mclapply(1:length(working.cell.split), Barcode_Calling_Function, mc.cores = numCores) #Match valid barcodes to cells
Barcodes.aligned.df<-data.frame(CBC=names(working.cell.split),VBC=unlist(Barcodes.aligned)) #Turn matched VBC/CBCs into a dataframe
Barcodes.aligned.df$VBC_Count<-as.numeric(lapply(Barcodes.aligned.df$CBC,Multi_VBC_count_function))
Barcodes.aligned.df$Tier<-ifelse(Barcodes.aligned.df$VBC==0,0,1)

#Making Final Barcode Index
Final.VBC.Index.df<-data.frame(VBC_Final=1:length(Full_VBC.Single.ls),Type=c(rep("MOI_1",length(MOI_1_VBC)),rep("SuperInfection",length(SuperInfection.VBCs.ls)),rep("SingleCellClone",length(VBC_Single_Cell.vec))))
Final.VBC.Index.df$Type<-as.character(Final.VBC.Index.df$Type)


#Getting Tier2 Barcodes
Multiplet.df<-Barcodes.aligned.df[Barcodes.aligned.df$VBC == 0,]


Dom.ls<-list()
for (i in unique(Multiplet.df$CBC)){
  working.set<-working.cell.split[[i]]
  barcode.set<-as.character(working.set$barcode)
  hits<-c()
  for (x in barcode.set){
    hits[x]<-grep(x,Full_VBC.Single.ls, fixed = TRUE)
  }
  combined.matches<-data.frame(VBC=hits)
  combined.matches$UMI_Count<-working.set$UMI_Count[working.set$barcode %in% rownames(combined.matches)]
  rownames(combined.matches)<-c()
  VBC.unique<-data.frame(VBC=unique(combined.matches$VBC))
  VBC.unique$UMI_Count<- lapply(VBC.unique$VBC, function(x) max(combined.matches$UMI_Count[combined.matches$VBC == x]))
  Dom.ls[[i]]<-VBC.unique
} #For Superinfection VBCs, will output max UMI count for highest expressed individual barcode

Tier2_exploration<-data.frame(CBC=names(Dom.ls))
Tier2_exploration$CBC<-as.character(Tier2_exploration$CBC)
Tier2.results.ls<-list()
for (i in Tier2_exploration$CBC){
  current.set<-Dom.ls[[i]]
  current.set$UMI_Count<-as.numeric(current.set$UMI_Count)
  Max.UMI<-max(current.set$UMI_Count)
  Sum.UMI<-sum(current.set$UMI_Count)
  Second.UMI<-as.numeric(current.set$UMI_Count)[order(as.numeric(current.set$UMI_Count))][length(current.set$UMI_Count)-1]
  Dom.bc<-current.set$VBC[current.set$UMI_Count==Max.UMI]
  results.merged<-paste(Max.UMI,Second.UMI,Sum.UMI,Dom.bc,sep="-")
  Tier2.results.ls[i]<-ifelse(length(results.merged)==1,results.merged, "0-0-0-0")
} #will output "0-0-0-0" in the event 2 barcodes occur with the exact same UMI
Tier2.results.df <- data.frame(do.call(rbind, str_split(Tier2.results.ls,"-")))
Tier2.results.df<-cbind(Tier2_exploration,Tier2.results.df)
names(Tier2.results.df)[2:5]<-c("Max","Second","Total","DOM_BC")
Tier2.results.df$Max<-as.numeric(as.character(Tier2.results.df$Max))
Tier2.results.df$Second<-as.numeric(as.character(Tier2.results.df$Second))
Tier2.results.df$Total<-as.numeric(as.character(Tier2.results.df$Total))
Tier2.results.df$RNA<-object.integrated$nCount_RNA[Tier2.results.df$CBC]

Tier2.results.df$Prop_Dom<-Tier2.results.df$Max/Tier2.results.df$Total
Tier2.results.df$Prop_Second<-Tier2.results.df$Max/Tier2.results.df$Second

Equal.UMI.set.df<-Tier2.results.df[Tier2.results.df$Max==0,] #create subset Tier2 results with equal UMI
Tier2.results.df<-Tier2.results.df[!Tier2.results.df$CBC %in% Equal.UMI.set.df$CBC,] #remove subset Tier2 results with equal UMI


Tier2.results.df$DOM_BC<-as.numeric(as.character(Tier2.results.df$DOM_BC)) #Convert DOM_BC from factor into interger for VBC type matching
Tier2.Viral.Subtype.vec<-c() #create empty vector in which to add VBC type in following line
for (i in 1:dim(Tier2.results.df)[1]){
  VBC.key<-Tier2.results.df[i,5]
  Output.type<-Final.VBC.Index.df[Final.VBC.Index.df$VBC_Final==VBC.key,2]
  Tier2.Viral.Subtype.vec[i]<-Output.type
} #add VBC type (ie. MOI_1, Superinfection, Singecellclone, etc.) for dominant VBC within in each cell
Tier2.results.df$VBC.type<-Tier2.Viral.Subtype.vec #add VBC type to Tier2 data.frame

ggplot(Tier2.results.df, aes(x=Prop_Second, y=Prop_Dom )) + geom_point() + theme_bw() + geom_vline(xintercept = 5, linetype="dotted", color = "red", size=1.5) +xlim(0,50) +ylim(0,1)
#The threshold for a Tier2 called cell can be set below. We used 5 for most of our analysis as the plot above seemed to indicate two groups of cells on either side of the cutoff. 
#For GW15_Rep2, we used a cutoff defined by the linear equation below which divided the data into two groups.
#Tier2.results.df$Tier2<-ifelse(Tier2.results.df$Prop_Dom > -0.05*Tier2.results.df$Prop_Second + 0.9,2,0) ## cutoff for Tier2 inclusion in GW15_Rep2 set



Tier2.results.df$Tier2<-ifelse(Tier2.results.df$Prop_Second >=5,2,0) #Standard cutoff for Tier2 inclusion (First/Second VBC UMI ratio >5)
Equal.UMI.set.df$VBC.type<-"Same_UMI_Count"
Equal.UMI.set.df$Tier2<-0 #Modify the set of potential Tier2 CBCs with equal counts of first/second VBCs in order to merge with rest of Tier2 DF
Tier2.results.df<-rbind(Tier2.results.df,Equal.UMI.set.df) #Merge Tier2 dataframe with Tier2 calls with subset of CBCs with equal #s of first and second highest VBC counts


Tier2.results.df$DOM_BC<-as.numeric(Tier2.results.df$DOM_BC)
Tier2.results.merge.df<-Tier2.results.df[Tier2.results.df$Tier2==2, c(1,5,9,10)]
colnames(Tier2.results.merge.df)<-c("CBC","VBC","VBC.type","UMI")

#Update Barcode DF with CBCs resolved as being Tier2
for (i in Barcodes.aligned.df$CBC[Barcodes.aligned.df$Tier==0]){
  position<-as.numeric(rownames(Barcodes.aligned.df[Barcodes.aligned.df$CBC== i, ]))
  Barcodes.aligned.df[position,2]<-ifelse(Tier2.results.df$Tier2[Tier2.results.df$CBC == i] ==2, Tier2.results.df$DOM_BC[Tier2.results.df$CBC == i], 0)
  Barcodes.aligned.df[position,4]<-ifelse(Tier2.results.df$Tier2[Tier2.results.df$CBC == i]==2, 2, 0)
}

#Add Multiplet entry to Final.VBC.Index.df
Final.VBC.Index.df<-rbind(c(0,"Multiplet"),Final.VBC.Index.df)
Final.VBC.Index.df$VBC_Final<-as.numeric(Final.VBC.Index.df$VBC_Final)


#Add VBC type to Final.VBC.Index.df
Viral.Subtype.vec<-c()
for (i in 1:dim(Barcodes.aligned.df)[1]){
  VBC.key<-Barcodes.aligned.df[i,2]
  Output.type<-Final.VBC.Index.df[Final.VBC.Index.df$VBC_Final==VBC.key,2]
  Viral.Subtype.vec[i]<-Output.type
}
Barcodes.aligned.df$VBC.type<-Viral.Subtype.vec

Barcodes.aligned.df$CBC<-as.character(Barcodes.aligned.df$CBC)

#Re-adjust Tier status of Superinfection clones so that cells that express all barcodes that define a clone are considered tier 1, and those that express a subset are tier 2
SuperInfection.subset.df<-Barcodes.aligned.df[Barcodes.aligned.df$VBC.type=="SuperInfection",]

SuperInfection.Tier.ls<-list()
for (i in unique(SuperInfection.subset.df$CBC)){
  working.set<-working.cell.split[[i]]
  barcode.set<-as.character(working.set$barcode)
  hits<-c()
  for (x in barcode.set){
    hits[x]<-grep(x,Full_VBC.Single.ls, fixed = TRUE)
  }
  combined.matches<-data.frame(VBC=hits)
  Final.subset<-SuperInfection.subset.df[SuperInfection.subset.df$CBC ==i,]
  VBC_called<-unlist(Full_VBC.Single.ls[Final.subset$VBC])
  Intersection.length<-length(intersect(VBC_called,rownames(combined.matches)))
  Total.VBC.length<-length(VBC_called)
  Final.SuperInfection.Tier<-ifelse(Intersection.length==Total.VBC.length,1,2)
  SuperInfection.Tier.ls[[i]]<-Final.SuperInfection.Tier
} 
SuperInfection.Tier.df<-data.frame(CBC=unique(SuperInfection.subset.df$CBC),Tier=unlist(SuperInfection.Tier.ls))

for (i in SuperInfection.Tier.df$CBC){
  Barcodes.aligned.df[Barcodes.aligned.df$CBC == i,4]<-SuperInfection.Tier.df$Tier[SuperInfection.Tier.df$CBC==i]
} #Update VBC/CBC list with SuperInfection Tiers


#Save Temp Files
saveRDS(Barcodes.aligned.df,"/outpath/Sample_BC_calls.rds")
saveRDS(Full_VBC.Single.ls,"/outpath/Sample_Full_VBC.Single.ls.rds")
saveRDS(Final.VBC.Index.df,"outpath/Sample_Final.VBC.Index.df.rds")

#Create vector of all samples
FirstRoundSet.set<-c("Sample1","Sample2","ETC")


#Make final STICR/10X Cell barcode data.frame
BC_path<-"/outpath/"
merged.bc.ls<-list()
for (i in 1:length(FirstRoundSet.set)){
  working.val<-FirstRoundSet.set[i]
  path1<-paste0(working.val,"_BC_calls.rds")
  path2<-paste0(working.val,"_Full_VBC.Single.ls.rds")
  path3<-paste0(working.val,"_Final.VBC.Index.df.rds")
  Barcodes.aligned.df<-readRDS(paste0(BC_path,path1))
  Full_VBC.Single.ls<-readRDS(paste0(BC_path,path2))
  Final.VBC.Index.df<-readRDS(paste0(BC_path,path3))
  for (j in 1:dim(Barcodes.aligned.df)[1]){
    working.set<-Barcodes.aligned.df[j,]
    updated.val<-ifelse(working.set$VBC == 0, "Multiplet", ifelse(working.set$VBC.type=="SuperInfection", paste(unlist(Full_VBC.Single.ls[working.set$VBC]), collapse =":"),Full_VBC.Single.ls[working.set$VBC]))
    Barcodes.aligned.df[j,6]<-updated.val
  }
  colnames(Barcodes.aligned.df)[6]<-"Sequence"
  Barcodes.aligned.df$Library<-FirstRoundSet.set[i]
  merged.bc.ls[[i]]<-Barcodes.aligned.df
}

merged.bc.df<-do.call(rbind.data.frame, merged.bc.ls)

#save Final Dataset
saveRDS(merged.bc.df,"/pathout/Final_VBC_CBC_merged.rds")


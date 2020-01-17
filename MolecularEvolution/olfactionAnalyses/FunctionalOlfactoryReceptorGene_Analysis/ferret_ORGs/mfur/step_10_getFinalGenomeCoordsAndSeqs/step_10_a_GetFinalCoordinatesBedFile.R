### Have ORF Info from Step 3d
### Have Starting M infrom from step 5d

### Have coordinates of hit from name 
### Bed file was +-300 of that. [If you change this you need to know]
spp="mfur"
setwd(paste("/Users/annabelbeichman/Documents/UCLA/Otters/Olfaction/20170707_MainOlfactionResults/",spp,sep=""))
# this is the original bed file that you pulled sequences out of (with coords +- 300bp)
bed <- read.table(paste("step_1_c_CleanedByEvalue_ABScript/",spp,".cleanedByEvalue.LengthFilterOnly.Pruned.scaff.0based.start.minus300.stop.plus300.name.eval.strand.bed",sep=""))
head(bed)
colnames(bed) <- c("Scaffold","BedStartMinus300","BedStopPlus300","bedID","evalue","strand") # these are zero based!
#### These are the results of picking ORFs in step 3
orfCoords <- read.table("step_3_findORFs/ORF_longThan_810_bp_step_2_results.olfacUniqueBlastHits.stranded.fasta.ORFStStopTable",header=F)
head(orfCoords)
colnames(orfCoords) <- c("bedID","frame","orfStart","orfStop") # these are 1 based
# need to make them 0 based (don't chang orfStop, just subtract 1 from orfStart to be bed compatible) (stop coordinate is non inclusive)
orfCoords_0Based <- orfCoords
orfCoords_0Based$orfStart <- orfCoords$orfStart - 1
head(orfCoords_0Based)
### These are the location of the best "M" AA within the **protein** sequence
# NOTE: genome coords (bed) and orfCoords are in terms of NUCLEOTIDES WITHIN THE SEQUENCE
# But the M coord is in terms of AMINO ACIDS!
MCoords <- read.table("step_5_chooseM/step_5_result.pickedM.CoordinatesOfBestStart.txt",header=F)
head(MCoords)
colnames(MCoords) <- c("bedIDPlusAA","MCoord")
MCoords$MCoord_nt <- MCoords$MCoord *3
### Try merging stuff
bedOrf <- (merge(bed,orfCoords_0Based,by="bedID"))
head(bedOrf)
# BEFORE merging in MCoords, check number of aas to make sure they'll match:
### Am doing frame right
bedOrf$AACount <- (bedOrf$orfStop - bedOrf$orfStart)/3
bedOrf$bedIDPlusAA <- paste(bedOrf$bedID,"-",bedOrf$AACount,"_aa",sep="")
head(bedOrf)
bedOrfM <- merge(bedOrf, MCoords,by="bedIDPlusAA") # ends up being dimensions of M coords
## THIS IS SUPERIOR because it gets rid of multiple hits with diff reading frames
# not necessary in the other species I think.
# because some don't have an ORF, some were visually excluded before M picking
### For positive strand, to get coordinates 
### First! Deal with MCoord and Orf Coords.
# Should be straightforward. (don't need to treat +- separately at this stage)
bedOrfM$orfStartCorrectM <- bedOrfM$orfStart + bedOrfM$MCoord_nt   # has to be in terms of nucleotides!!!!! 

# now need to get actual genomic coordinates (0 based)
# going to be different for STRANDS
bedOrfM$finalBedStart <- NA
bedOrfM$finalBedEnd <- NA
# for positive strand: do bedStart + orfStartCorrectM to get final start position
# and bed start + orfEnd for final end position 
# for negative strand: do bedEnd - orfSt and bedend - orfEnd.
bedOrfM[bedOrfM$strand=="+",]$finalBedStart <- bedOrfM[bedOrfM$strand=="+",]$BedStartMinus300 + bedOrfM[bedOrfM$strand=="+",]$orfStartCorrectM
bedOrfM[bedOrfM$strand=="+",]$finalBedEnd <- bedOrfM[bedOrfM$strand=="+",]$BedStartMinus300 + bedOrfM[bedOrfM$strand=="+",]$orfStop # I think this doesn't change based on M because it's all relative to the original sequence pulled out of the genome (taht is between bedstart and bedend)
head(bedOrfM)

bedOrfM[bedOrfM$strand=="-",]$finalBedEnd <- bedOrfM[bedOrfM$strand=="-",]$BedStopPlus300 - bedOrfM[bedOrfM$strand=="-",]$orfStart
bedOrfM[bedOrfM$strand=="-",]$finalBedStart <- bedOrfM[bedOrfM$strand=="-",]$BedStopPlus300 - bedOrfM[bedOrfM$strand=="-",]$orfStop

# also want to get final bed IDs (and somehow log original IDs?)
finalBed <- bedOrfM[,c("Scaffold","finalBedStart","finalBedEnd","bedIDPlusAA","evalue","strand")]
head(finalBed)

write.table(finalBed,"step_10_getFinalGenomeCoordsAndSeqs/finalORCoordinates.all.0based.bed",quote=F,sep="\t",row.names = F,col.names = F)

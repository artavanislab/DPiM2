module add R

SCRDIR="/home/glocke/DPiM/scripts"

INPUT=$1
OUTPUT=$2

DATA=$INPUT.cp

# make a backup copy
cp $1 $DATA

# replace <97> special characters
perl -pi -e 's/\x97/-/g' $DATA

# remove 'FH' from the Tap_ids such as FHEGFP -> EGFP
perl -pi -e 's/FH//g' $DATA


# (total number of ## protein ids)/(total protein ids - ## proteins ids) to estimate protein FDR
#Rscript $SCRDIR/checkFDR.R $DATA

# remove lines contain '##FBgn' (ids such as ##FBgnxxx which are peptide matches
# to the reverse translation of the genome (called decoys) and it is a way
# Gygi lab calculates the false discovery rate of peptide matches in the mass
# spec data. It is safe it ignore them for the HGSCore analysis.

sed -i '/##FBgn/d' $DATA

# remove all rows with "_TAG" 
# The "_TAG" corresponds to the peptide with FLAG-HA sequence and some bit 
# of the bait protein sequence. We usually add it that to the corresponding 
# FBGn, which is safe to do for the HGSCore analysis. Julian added:
# these _TAG ids should be removed entirely because they are unreliable

#perl -pi -e 's/_TAG//g' $DATA
sed -i '/_TAG/d' $DATA

#sed -i 's/Bob Obar/Bob_Obar/' $DATA
#sed -i 's/Pujita Vaidya/Pujita_Vaidya/' $DATA
#sed -i 's/Christina Wong/Christina_Wong/' $DATA

# replace 'FBgn0266084' by 'FBgn0261259'
# FBgn0266084 corresponds to gene called Fhos, which also have other secondary FBgn IDs listed below:
# FBgn0010755
# FBgn0035925
# FBgn0035927
# FBgn0040229
# FBgn0052025
# FBgn0052030
# FBgn0261259
# only FBgn0261259 is found in our dataset, then replace FBgn0266084 with that

perl -pi -e 's/FBgn0266084/FBgn0261259/g' $DATA

# change the file from dos to unix format
#dos2unix $DATA

#Rscript
#Rscript $SCRDIR/parseDPiMdata.r $DATA $OUTPUT
mv $DATA $OUTPUT

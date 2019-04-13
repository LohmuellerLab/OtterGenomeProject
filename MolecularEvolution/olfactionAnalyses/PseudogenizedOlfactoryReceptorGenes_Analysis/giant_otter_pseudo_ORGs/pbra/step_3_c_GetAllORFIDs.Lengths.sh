# get names:

spp=$1
result=$2

# fasta_tool adds extra space to the name, get rid of it

sed -i'' 's/ //g' $result # because it went through fasta_tool it got an extra space
# get rid of "found gene"

sed -ibak -E 's/[0-9]+genefound//g' $result
# pull out all ORFs (above whatever min threshold) 
grep ">${spp}" $result > ${result%.fasta}.IDlist

sed -i'' 's/ //g' ${result%.fasta}.IDlist  # because some results will and won't have this depending on 
# if they went through fasta_tool. this makes sure it's uniform

# then put it into R 

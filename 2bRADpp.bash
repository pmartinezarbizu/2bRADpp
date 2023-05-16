#!/bin/bash  
### Pipeliine for 2bRAD processing of raw reads
###use bbmap from https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/

# cd ~/projects/RADs/PipelineRads/
# run like: bash ./2bRADpp.bash  > logfile.txt 2>&1 &

#####################################
### EDIT THESE OPTIONS            ###
#####################################

	# your bbmap path
	bbpath='/home/pmartinez/Downloads/bbmap'

	# define enzyme pattern uncomment 4 lines per enzyme
	
	# BcgI
	EZN=BcgI
	ENZ1U=.G.{10}CGA.{6}TGC.{10}C.
	ENZ1L=.G.{10}GCA.{6}TCG.{10}C.	
	ENZ2U=.{10}CGA.{6}TGC.{10}
	ENZ2L=.{10}GCA.{6}TCG.{10}
 	
	
	# BsaXI	
	#EZN=BsaXI	
	#ENZ1U=..G.{9}AC.{5}CTCC.{8}C. 	
	#ENZ1L=.G.{8}GGAG.{5}GT.{9}C.. 
	#ENZ2U=.{12}AC.{5}CTCC.{10} 
	#ENZ2L=.{10}GGAG.{5}GT.{12} 	
	

#	Adapter1	GTCT
#	Adapter2	AGAC
#	Adapter3	ACCA
#	Adapter4	AGTG
#	Adapter5	CATC

	# define adapter used
	#BC=AGTG

	#declare an array of adapter barcode sequences
	declare -a bcode=('GTCT' 'AGAC' 'ACCA' 'AGTG' 'CATC')
	#declare -a bcode=('CATC')


# check format must be: 1s1_L001_R1_.fastq

#L001 is third position see line 64:  l=$(cut -d_ -f2 <<< "$f")

#####################################
### DO NOT EDIT BEYOND THIS LINE  ###
#####################################
	
# start
	#echo date "Logfile batch RADs processing" > $EZN.$BC.logfileRads.txt
	#echo Enzyme is $EZN, pattern $ENZ1U >> $EZN.$BC.logfileRads.txt
	#echo Barcode of Adapter ist $BC >> $EZN.$BC.logfileRads.txt



for f in *_R1_*.fastq.gz; do
	    r=$(sed -e "s/_R1_/_R2_/" <<< "$f")
	    s=$(cut -d_ -f1 <<< "$f")
	    l=$(cut -d_ -f3 <<< "$f")	

	echo
	echo ====================================
	echo Processing sample $s
	echo ====================================


	echo ====================================	
	echo "Step 1 : discover adapters, use paired-end reads"
	echo ====================================
	echo
	$bbpath/bbmerge.sh in1=$f in2=$r outa=adapters.fa
	echo cat adapters.fa
	echo

	echo ====================================
	echo "Step 2: trim adapters and make PE contigs, filter sequences with average Q < 30 "
	echo ====================================

	#$bbpath/bbmerge.sh in1=$f in2=$r out=merged.$s.$l.fq outu=unmerged.$s.$l.fq \
	#ihist=ihist.txt adapters=adapters.fa minavgquality=10
	
	$bbpath/bbmerge.sh in1=$f in2=$r \
	out=merged.$s.$l.fq outu=unmerged.$s.$l.fq \
	rem qout=33 maxcalledquality=39 extend2=44 k=16 ecco mix

	echo ====================================
	echo Step 3: remove duplicate sequences, considering also reverse complements.
	echo Here the first 4 bases detecting PCR duplicates are still retained,
	echo But PCR duplicates are removed.
	echo "Also the last 4 barcode-bases are retained"
	echo ==================================== 
	

	$bbpath/dedupe.sh in=merged.$s.$l.fq out=mergeddedup.$s.$l.fq overwrite=true
	#$bbpath/dedupe.sh in=$f out=mergeddedup.$s.$l.fq overwrite=true

	echo ==================================== 
	echo "
	Step 4: demultiplex by barcode,  keep sequence name....\n
	 	"
	echo ==================================== 

	for BC in "${bcode[@]}"; do
	
	echo Processing Enzyme $EZN, pattern $ENZ1U with Adapter $BC	

	#extract upper strand	
	grep -P -B1 --no-group-separator $ENZ1U$BC mergeddedup.$s.$l.fq > $s.$l.up.fq
	
	#extract lower strand 
	grep -P -B1 --no-group-separator $ENZ1L$BC mergeddedup.$s.$l.fq > $s.$l.low.fq


	# Step 5: Extract sequence names
	grep -P --no-group-separator '^@' $s.$l.up.fq > $s.$l.up.names.fq
	grep -P  --no-group-separator '^@' $s.$l.low.fq > $s.$l.low.names.fq
	cup=$(grep -Pc '^@' $s.$l.up.fq)
	clow=$(grep -Pc '^@' $s.$l.low.fq)
	
	
	# Step 6: cut the exact RAD excluding barcodes 
	#grep -Po '.{12}CGA.{6}TGC.{12}' $s.up.fq > $s.up.seq.fq
	#grep -Po '.{12}GCA.{6}TCG.{12}' $s.low.fq | tr GATCgatc CTAGCTAG | rev > $s.low.seq.fq
	grep -Po --no-group-separator $ENZ2U $s.$l.up.fq > $s.$l.up.seq.fq
	grep -Po --no-group-separator $ENZ2L $s.$l.low.fq | tr GATCgatc CTAGCTAG | rev > $s.$l.low.seq.fq

	#Step 7: Paste Names and Rads into one file
	paste -d '\n' $s.$l.up.names.fq $s.$l.up.seq.fq > $s.$l.up.final.fq
	paste -d '\n' $s.$l.low.names.fq $s.$l.low.seq.fq > $s.$l.low.final.fq
	
	#Step 8: convert final result to fasta format
	cat $s.$l.up.final.fq $s.$l.low.final.fq  | sed 's/^@/>/' > $s.$l.$EZN.$BC.fasta
	cfinal=$(grep -Pc '^>' $s.$l.$EZN.$BC.fasta)
	#cfinalgood=$(grep -Pc '.{12}CGA.{6}TGC.{12}' $s.final.fasta)
	cfinalgood=$(grep -Pc $ENZ2U $s.$l.$EZN.$BC.fasta)
	
	echo == 
	echo Sample $s.$l. has $cup RADs on upper strand and $clow RADs on lower strand 
	echo  
	echo Sample $s.$l. has $cfinal RADs in total
	echo Sample $s.$l. has $cfinalgood good RADs 
	echo == 	
	echo Sample $s.$l. has $cfinalgood good RADs >> $EZN.$BC.logfileRads.txt 
	
# copy to final file
	cat $s.$l.$EZN.$BC.fasta >> $s.$EZN.$BC.fasta
	rm $s.$l.$EZN.$BC.fasta
	
	done #done barcodes

	#Step 10: cleaning up ...
	rm merged.$s.$l.fq

	rm $s.$l.up.fq
	rm $s.$l.low.fq

	rm $s.$l.up.names.fq
	rm $s.$l.low.names.fq

	rm $s.$l.up.seq.fq
	rm $s.$l.low.seq.fq

	rm $s.$l.up.final.fq
	rm $s.$l.low.final.fq

	rm mergeddedup.$s.$l.fq
done
	#delete empty files	
	find . -type f -size 0 -delete
	echo ...done




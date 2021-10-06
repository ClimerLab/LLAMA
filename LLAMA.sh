#!/bin/bash

GML_FILE=/media/ken/ExtraSpace/School/Research/DataSets/networks/football.gml
RESULTS_DIR=/media/ken/ExtraSpace/School/Research/LLAMA/
OUTPUT_TAG=football_steadystate_Q
NUM_ISLANDS=30
CFG_FILE=LLAMA_Q.cfg

OUTPUT_DIR=$RESULTS_DIR$OUTPUT_TAG/
STATS_DIR="$OUTPUT_DIR"Stats/

#---------------------------------------------------------------------------------------------------
# Seperate graph into components
#---------------------------------------------------------------------------------------------------
if [ -f "NUM_COMP.txt" ]; then
	rm "NUM_COMP.txt"
fi

# Check for dir, if not found create it using mkdir #
if [ ! -d $OUTPUT_DIR ]; then
	mkdir -p $OUTPUT_DIR
fi

# Check for dir, if not found create it using mkdir #
if [ ! -d "$OUTPUT_DIR"Best/ ]; then
	mkdir -p "$OUTPUT_DIR"Best/
fi

./SepComp $GML_FILE $OUTPUT_DIR $OUTPUT_TAG
if [ -f $OUTPUT_DIR"NUM_COMP.txt" ]; then
	source "$OUTPUT_DIR"NUM_COMP.txt
	rm $OUTPUT_DIR"NUM_COMP.txt"
else
	NUM_COMP=-1
	TOTAL_NODES=0
	TOTAL_EDGES=0
fi

if [ $NUM_COMP -eq -1 ]; then
	echo "Error seperating graph...exiting"
	exit -1
fi

#---------------------------------------------------------------------------------------------------
# Run LLAMA
#---------------------------------------------------------------------------------------------------

# Check for dir, if not found create it using mkdir #
if [ ! -d $STATS_DIR ]; then
	mkdir -p $STATS_DIR
fi

for ((num=1; num<=$NUM_ISLANDS; num++))
do
	if [ ! -d "$OUTPUT_DIR"I_$num ]; then
		mkdir -p "$OUTPUT_DIR"I_$num
	fi	
done


# Loop through all compoents found by "SepComp"
for ((CUR_COMP=1; CUR_COMP<=$NUM_COMP; CUR_COMP++))
do
	echo "Processing comp $CUR_COMP"
	# Reset the number of completed islands for current component
	COMPLETED_ISLANDS=0
	# Each island will iterate through a given number of generations and then a migration event will occur. This while loop
	# calls LLAMA on each island until all generations have been completed
	while [ $COMPLETED_ISLANDS -lt $NUM_ISLANDS ]
	do
		# Loop through all islands for the current set of iterations
		for ((CUR_ISLAND=1; CUR_ISLAND<=$NUM_ISLANDS; CUR_ISLAND++))
		do
			# Check if GML files exists
			if [ -f "$OUTPUT_DIR$OUTPUT_TAG"_comp"$CUR_COMP".gml ]; then
				# Check if clustering result file exists
				if [ -f "$OUTPUT_DIR"I_"$CUR_ISLAND/$OUTPUT_TAG"_comp"$CUR_COMP"_"$CUR_ISLAND".out ]; then
					# Increment the number of completed islands for this component
					((COMPLETED_ISLANDS = COMPLETED_ISLANDS + 1))
				else
					# Call LLAMA
					./LLAMA "$OUTPUT_DIR$OUTPUT_TAG"_comp"$CUR_COMP".gml $CFG_FILE $TOTAL_NODES $TOTAL_EDGES $OUTPUT_DIR $OUTPUT_TAG $CUR_ISLAND $STATS_DIR #$NUM_ISLANDS
				fi
			fi
		done	
	done

	# Delete tmp files from folders
	for ((CUR_ISLAND=1; CUR_ISLAND<=$NUM_ISLANDS; CUR_ISLAND++))
	do
		if [ -f "$OUTPUT_DIR"I_"$CUR_ISLAND"/popfile.txt ]; then
			rm "$OUTPUT_DIR"I_"$CUR_ISLAND"/popfile.txt
		fi
		if [ -f "$OUTPUT_DIR"I_"$CUR_ISLAND/$OUTPUT_TAG"_init_popfile.txt ]; then
			rm "$OUTPUT_DIR"I_"$CUR_ISLAND/$OUTPUT_TAG"_init_popfile.txt
		fi	
	done
done


#---------------------------------------------------------------------------------------------------
# Find best cluster assignment for each component
#---------------------------------------------------------------------------------------------------
./FindBestComp $OUTPUT_DIR $OUTPUT_TAG $NUM_COMP $NUM_ISLANDS


#---------------------------------------------------------------------------------------------------
# Combine best cluster assignments
#---------------------------------------------------------------------------------------------------
SINGLETON_FILE="$OUTPUT_DIR$OUTPUT_TAG"_comp0.nn
./CombComp "$OUTPUT_DIR"Best/ $OUTPUT_TAG $GML_FILE $SINGLETON_FILE $OUTPUT_DIR $NUM_COMP

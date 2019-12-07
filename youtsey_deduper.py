#!/usr/bin/env python3
import os
import argparse
import re

########MANAGE INPUTS
def get_args():
    parser = argparse.ArgumentParser(description="Deduplicates SAM File of single-end alignments. Outputs the input SAM file name with (_deduped.sam) extension")
    parser.add_argument("-f", "--file", help= "SAM file to be deduplicated (REQUIRED)", required = True)
    parser.add_argument("-p", "--paired", help= "Specify if paired-end reads aligned (DO NOT USE: NOT DEVELOPED YET)", required = False, default= False)
    parser.add_argument("-u", "--umi", help= "File containing column of UMI's. If no UMI file is provided, the script will consider every unique UMI encountered (OPTIONAL)", required = False)


    return parser.parse_args()

args = get_args()
inputSAM = args.file
umiFile = args.umi

#This script has two main parts.
#First (In Python): Find and calculate the true left most position (accounting for: soft-clipping, insertions, deletions, gaps, & strandedness)
                  # concat _UMI+trueLeftMost_ to the end of the first field for each read
                  # save file to remove duplicates of _UMI+trueLeftMost+strand_ (save header separately)
#Second (In os.system): run sort command on _UMI+trueLeftMost_ and keep only unique ID's encountered first
                      # concatenate header back onto deduplicated sam file
                      # remove temporary files (undeduped SAM body and SAM header files)
###########READ IN UMI's from UMI FILE

#tuple to store all UMI's. Should be left empty if no umi file. If empty, the SAM line will be processed regardless of UMI sequence
umiList = list()
if umiFile is not None:

    with open(umiFile, "r") as umiInput:
        for line in umiInput:
            umiList.append(line.strip())

###########DEFINE FUNCTIONS
def getUMI(header_string):
    '''Inputs first field in SAM row. Returns the UMI in the last : field'''
    fields = header_string.split(":")
    return fields[len(fields)-1]

def createCustomFlag(umi,leftmost, strand):
    '''Inputs umi, left-most position, & strandedness. Returns _UMI+trueLeftMost+strand_ string'''
    return "_"+umi+ "+" + str(leftmost)+ "+" + strand+ "+" + strand + "_"



###########RUN MAIN PROGRAM BODY
#PART ONE (Python)
with open(inputSAM,"r") as samFile, open("tmp.header.sam", "w") as headerOut, open("tmp.sam", "w") as samOut:


    for line in samFile:

            #check if header line:
            if line[0] == "@":
                #Save to header file. Will be combined back with tmp.sam after sort -u
                headerOut.write(line)
            else:
            #Not header line. Calculate trueLeftMost and append _UMI+trueLeftMost+strand_
                flags = line.split("\t")
                #will be combined with the reformatted header later
                afterHeader = re.findall(r"(\t.*)", line)[0]

                header = flags[0]
                umi = getUMI(header)
                bitFlag = flags[1]
                cigarString = flags[5]
                #original left-most position (not accounting for soft clipping)
                orig_leftmost = int(flags[3])
                adjusted_leftmost = orig_leftmost

                #checking to see if UMI's match list or UMI's were not specified
                if umi in umiList or len(umiList) == 0:

                    #adjusting for forward read
                    if cigarString[1] == "S":
                        adjusted_leftmost = adjusted_leftmost - int(cigarString[0])

                    customFlag = createCustomFlag(umi, adjusted_leftmost, "F")
                    samOut.write(header+ customFlag + afterHeader +"\n" )





#PART TWO (os.system)
#sort command will provide new SAM file with only unique _UMI+trueLeftMost+strand_ (based on first encountered)
commandString = "sort -t _ -k2,2 -u tmp.sam > "

#outputFile will encorporate the input file name
inputFileName = inputSAM.split(".")[0]
outputFileString = inputFileName + "_deduped.sam"

#save deduped SAM body to temp file
os.system(commandString + "tmp.deduped.sam")

#recombine SAM body file with header file
os.system("cat "+ "tmp.header.sam tmp.deduped.sam > "+ outputFileString )

#remove temporary files
os.system("rm tmp.deduped.sam tmp.header.sam tmp.sam")

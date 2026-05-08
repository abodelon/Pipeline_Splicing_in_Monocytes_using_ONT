#!/usr/bin/env python

#############################
# FLAIR database generation #
#############################

import os, sys, string, HTSeq, time

# Arguments
dataDirBase= sys.argv[1] 
dataDirBaseBulk= sys.argv[2] 
folderName= sys.argv[3] 
folderNameSamples= sys.argv[4] 
folderNameBatch= sys.argv[5] 
outputFN= sys.argv[6] 

def makeMonocytesOnlyNGSbasedJunctionsFromSamFile(dataDirBase=None, dataDirBaseBulk=None, folderName=None, folderNameSamples=None, folderNameBatch=None, outputFN=outputFN):
	# the junctions_from_sam.py program from FLAIR, is unable to capture the strand specific information from the STAR alignment files
	# Therefore, I am getting the junctional information out myself
	# pipeline: 
	# - go through each alignement
	# - select mono-aligned reads (uniquely aligned), with a deletion (from genome, to get to read) > 100 base-pairs, or an N-part in the cigar
	#   . Also, deletion/N part must be flanked by match positions at least 5 base-pairs long on both sides
	# - store begin/end/strand/chrom info of the N-part or of the deletion (introns)
	# - keep introns that are found in at least 2 samples
	# - return kept introns to bed
	# Another problem 1:
	# - lots of reads are marked inconsistent because an exon is not annotated and not found in the NGS (eg gene LGALS3)
	#   . The problem still only exists with unannotated splice donor/acceptor sites, even if the current junction is not annotated it can be
	#     consistent by 2 other junctions (ie skipped exon is allowed)
	# - Possible solution 1: use mapped ONT data to obtain junctions shared by at least 2 samples (same as for NGS)
	#   . pro: all relevant junctions are found
	#   . con: junctions are based on the data that will later be analyzed
	#   . Ask for support of at least 1 NGS sample, and 2 nanopore sample donors (IMPLEMENTED)
	#   . Ask for support from the nanopore samples of all (3) donors (IMPLEMENTED)
	# Another problem 2: 
	# - small exons that are mis-aligned to the begin/end of the exon after/before it (eg gene MAPK1IP1L) 
	# - Possible solution : Restore alignment to better matching alternative
	#   . Seek for mis-aligned overhangs of an exon, with a possible alternative exon before it
	#   . Look for much better alignment, if so restore
	#   . pro: get better alignments
	#   . con: too much work to implement 
	# - See function correctingMinimap2alignedDirectRNAseqBam (not IMPLEMENTED, too much work/too difficult)
	# Another problem 3: 
	# - exon/exon junctions that are not annotated and not present in NGS data (eg gene SOS2)
	#   . possible solution same as for the problem 1 solution
	
	if not dataDirBase:
		raise Exception('add where is your ONT alignment path')

	if not dataDirBaseBulk:
		raise Exception('add where is your NGS alignment path')

	# get the metadata from the NGS alignment
	perDatasetMetadata = giveAllPerDatasetMetadata_sysJIAstudy(folderName=folderName, folderNameSamples=folderNameSamples, folderNameBatch=folderNameBatch)
	# only get the information of the monocytes
	allMonoHCsamples = [smp for smp in perDatasetMetadata.keys() if (perDatasetMetadata[smp]['tissue'] == 'mono' and
																	 perDatasetMetadata[smp]['disease'] == 'healthy' and
																	 perDatasetMetadata[smp]['treatment'] in ['exvivo',
																											  'LPS'])]
	# get the name of the ONT samples
	allMonoFiles=os.listdir(dataDirBase)
	# Select the specific bam filles
	allMonoONTsamples = [x for x in allMonoFiles if (".bam" in x) and (not ".bai" in x)]

	StrandSwitcher = {'+': '-', '-': '+'}

	# create the new variables
	perIntronSamples = {}
	perIntronSamples_ONTonly = {}
	# save the ids of each nanopor long sequencing
	ONTdonorNs = {'PIC3822': 'donorM1', 'PIC3823': 'donorM1', 'PIC4991': 'donorM2', 'PIC4992': 'donorM2',
				  'CAL8276': 'donorM3', 'CAL8302': 'donorM3'}
	for NGSorONT, theSamples in [('NGS', allMonoHCsamples),
								 ('ONT', allMonoONTsamples)]:  # enter to all samples (ONT and NGS)
		for smp in theSamples:  # for each sample
			if NGSorONT == 'NGS':  # if it is NGS get the information from
				BAMfilename = dataDirBaseBulk + '/%sAligned.sortedByCoord.out.bam' % (
					smp)  # select the BAMfile
			if NGSorONT == 'ONT':  # if it is ONT
				# get the id of the sample (not the description)
				PICorCALn = \
					[j for i in smp.split('/') for j in i.split('_') if ('PIC' in j or 'CAL' in j) and (not 's' in j)][
						0]
				ONTdonorN = ONTdonorNs[PICorCALn.replace('.bam', '')]  # access to the donor id of the id of the sample
				BAMfilename = dataDirBase + "/" + smp  # prepare the bam file name
			if not os.path.isfile(BAMfilename):  # if the file does not exist, error
				raise Exception('fail')

			bam_reader = HTSeq.BAM_Reader(BAMfilename)  # read the bam file
			cnter = 0  # counter to 0
			for alignment in bam_reader:  # read each alinment line by line
				cnter += 1  # count
				if cnter % 100000 == 0:  # print the counter
					print(smp, " - sample", (allMonoONTsamples + allMonoHCsamples).index(smp), "from ",
						  len(allMonoONTsamples + allMonoHCsamples), ": ", cnter, "", '\r')
					sys.stdout.flush()  # print it one by one

				# check the alignments:
				if NGSorONT == 'NGS' and (not dict(alignment.optional_fields)['NH'] == 1):
					continue  # if it is not a uniquely aligned read in NGS, go to the next line
				if NGSorONT == 'ONT' and alignment.not_primary_alignment:
					continue  # if it is not a primary alignment, go to the next line

				# The 'CIGAR' (Compact Idiosyncratic Gapped Alignment Report) string is how the SAM/BAM format represents spliced alignments
				# Length cigar has to be more tha three: match - gap - match
				if len(alignment.cigar) >= 3:  # if the length of the cigar is mor ethan 3
					cigarPieces = alignment.cigar  # save the cigar
					for cigIndex in range(1, len(cigarPieces) - 1):  # enter each cigar
						if cigarPieces[cigIndex].type in ['N',
														  'D']:  # if there is a cigar type: ('D', 'deleted'), ('N', 'skipped')
							# if there is a gap of less than 100 bp: select introns of at least 100 bp
							if cigarPieces[cigIndex].type == 'D' and cigarPieces[cigIndex].ref_iv.length < 100:
								continue  # next loop
							# if before and after the N|D there is a match ('M', 'matched')
							if cigarPieces[cigIndex - 1].type == cigarPieces[cigIndex + 1].type == 'M':
								# and if the length of this match is at least 5 bp
								if cigarPieces[cigIndex - 1].ref_iv.length >= 5 and cigarPieces[
									cigIndex + 1].ref_iv.length >= 5:
									chrom = cigarPieces[cigIndex].ref_iv.chrom  # save the chromosome
									start = cigarPieces[cigIndex].ref_iv.start  # save the start
									end = cigarPieces[cigIndex].ref_iv.end  # save the end
									strand = cigarPieces[cigIndex].ref_iv.strand  # save the strand

									if NGSorONT == 'NGS':  # if it was NGS
										strand = StrandSwitcher[strand]  # switch the strand # complementary strand (?)

									intronInfo = (chrom, start, end, strand)  # and also the information of the intron

									# to check if the intron was found in two or more NGS 
									if NGSorONT == 'NGS':  # if it is NGS
										if not intronInfo in perIntronSamples:  # if it was not saved
											perIntronSamples[intronInfo] = set(
												[smp])  # add in the key of the info of the intron, the sample
										# if it was only found for now in one sample and the sample is different than the previous one
										elif len(perIntronSamples[intronInfo]) < 2 and (
												not smp in perIntronSamples[intronInfo]):
											perIntronSamples[intronInfo].add(smp)  # add the ONT sample

									# to check if the intron was found in all ONT donors
									if NGSorONT == 'ONT':  # if it is ONT
										if not intronInfo in perIntronSamples_ONTonly:  # and it is not described in perIntronSamples_ONTonly
											perIntronSamples_ONTonly[intronInfo] = set(
												[ONTdonorN])  # add the intron as a key and the sample as a value
										else:  # else
											perIntronSamples_ONTonly[intronInfo].add(
												ONTdonorN)  # add the ONT sample as a value

	# get all the unique/different donnors (there are 6 samples and they come from 3 different donnors)
	allONTdonors = set(ONTdonorNs.values())
	OUTfile = open(outputFN, 'w')  # wirte the output
	junctionsDone = set([])  # create the new variable to save the junctions

	# enter the information of each intron
	for intronInfo in perIntronSamples:
		if len(perIntronSamples[intronInfo]) >= 2:  # if the intron was found in at least 2 NGS samples
			chrom, start, end, strand = intronInfo  # get the information
			OUTfile.write(
				'\t'.join([chrom, str(start), str(end), 'noName', '0', strand, '\n']))  # write it in the output
			junctionsDone.add(intronInfo)  # and add the information in junctionsDone
		elif len(perIntronSamples[intronInfo]) == 1:  # if it was detected only in one NGS sample
			# Also, adding those junctions that are found in 1 NGS sample but 2 ONT samples
			#   . Ask for support of at least 1 NGS sample, and 2 nanopore samples (TO IMPLEMENT)
			# if it was also described in two donnors using ONT 
			if intronInfo in perIntronSamples_ONTonly and len(perIntronSamples_ONTonly[intronInfo]) >= 2:
				chrom, start, end, strand = intronInfo
				OUTfile.write('\t'.join([chrom, str(start), str(end), 'noName', '0', strand, '\n']))
				junctionsDone.add(intronInfo)  # add the information
	# for each intron that was described only in ONT and was not already added in the output (-junctionsDone)
	for intronInfo in set(perIntronSamples_ONTonly.keys()) - junctionsDone:
		if perIntronSamples_ONTonly[intronInfo] == allONTdonors:  # if all three donnors are included
			chrom, start, end, strand = intronInfo  # save the intron information
			OUTfile.write(
				'\t'.join([chrom, str(start), str(end), 'noName', '0', strand, '\n']))  # write it in the output
	OUTfile.close()  # close the connexion

def giveAllPerDatasetMetadata_sysJIAstudy(excludeSamples=None, folderName=None, folderNameSamples=None, folderNameBatch=None):
	# collect meta data per data set from the file ReferenceGenome/metadata_mono.txt
	# - tissue: CD4 / CD8 / mono / neutro 
	# - treatment: exvivo/LPS/medium
	# - disease: systemicJIA_active / systemicJIA_remission / polyJIA_active / healthy 
	# - patientID: ... tbd
	# - batch: ... tbd (possibly sequencing batch)
	# - Other meta data factors?
	#   . Sex: male/female
	# - treatment: CREBPP, DMSO, iBET, med, NAM9
	# - disease: HC (HD) or JIA (SF)
	# - donor (SF1-SF5, HD1-HD5)
	# - batch: from perStudySamples.. 

	# if not folder given raise an error
	if not folderName:
		raise Exception('add where is your metadataJIArnaSeq file')

	if not folderNameSamples:
		raise Exception('add where is your PerStudySamples file')

	if not folderNameBatch:
		raise Exception('add where is your PerStudyBatch file')

	# open the metadataJIArnaSeq.txt
	afile = open(folderName, 'r')  # open the metadata
	# save the first line of the file as headerLineCells, separated by tab and removing space and new line
	headerLineCells = afile.readline().replace('\n', '').replace(' ', '').split('\t')
	# if the header des not include sample name, disease, tissue, treatment and patient ID
	if not headerLineCells == ['sampleName', 'disease', 'tissue', 'treatment', 'patientID']:
		raise Exception('wrong file or layoutOfFile')  # raise an error

	# open the study-sample information
	afilesamples = open(folderNameSamples, 'r')
	perStudySamples={}  # create the dictionary

	while 1:  # read the file and create a dictionary with ID: study ID and as value, all samples.
		aline = afilesamples.readline()  # read a line
		if not aline:  # if it is empty
			afilesamples.close()  # close the connexion
			break  # out of the loop
		perStudySamples[aline.split('\t')[0]]=aline.split('\t')[1].replace('\n', '').split(',')  # create the dictionary

	perSampleNameRunID = {}  # create a new variable to save the information
	for runID in perStudySamples:  # for loop to each sample in perStudySample usign the runID (sample ID)
		for sampleName in perStudySamples[runID]:
			if sampleName in perSampleNameRunID:  # if the sample ID in
				raise Exception("double entry?")  # raise an exception because there is a douuble sample
			perSampleNameRunID[sampleName] = runID  # save the run id

	# open the study-run information
	afilesamples = open(folderNameBatch, 'r')
	perRunIDrunBatch={}  # create the dictionary

	while 1:  # read the file and create a dictionary with the run id as a key and batch id as a value
		aline = afilesamples.readline() # read a line
		if not aline: # if it is empty
			afilesamples.close() # close the connexion
			break # out of the loop
		perRunIDrunBatch[aline.split('\t')[0]]=aline.split('\t')[1].replace('\n', '')  # create the dictionary

	OUT = {}  # save the output
	while 1:  # loop to enter each loop
		aline = afile.readline()  # read each line
		if not aline:  # if there are no more lines
			afile.close()  # close the connexion and
			break  # go out of the loop
		# save in cell the line separated by tab and removing space and new line
		cells = aline.replace('\n', '').replace(' ', '').split('\t')
		# save sample, disease, tissue, treatment and patient id in cells
		sampleName, disease, tissue, treatment, patientID = cells
		runID = perSampleNameRunID[sampleName]  # enter the run id
		runBatch = perRunIDrunBatch[runID]  # enter the run batch id
		# give the patient id from the sample name
		patientID = givePatientIDfromSampleName(sampleName)

		# check if all samples are okay
		# if there is a disease that is not these ones
		if not disease in ['systemicJIA_active', 'systemicJIA_remission', 'polyJIA_active', 'healthy']:
			raise Exception('fail')  # raise an error
		# if there is a tissue that is not these ones
		if not tissue in ['CD4', 'CD8', 'mono', 'neutro']:
			raise Exception('fail')  # raise an error
		# if there is a treatment that is not these ones
		if not treatment in ['exvivo', 'LPS', 'medium']:
			raise Exception('fail')  # raise an error
		# if there is already a sample with this name (already saved in the OUT)
		if sampleName in OUT:
			raise Exception('double entry?')  # print an error
		# if it is a sample that is set to exclude (so excludeSamples is also not empty)
		if excludeSamples != None and sampleName in excludeSamples:
			continue  # continue to the next loop

		# add the sample information in the output: the name of the sample, the disease, tissue, treatment, the run, and the run batch
		OUT[sampleName] = {'sampleName': sampleName, 'disease': disease,
						   'tissue': tissue, 'treatment': treatment, 'patientID': patientID,
						   'runID': runID, 'runBatch': runBatch}
	return OUT  # return the output

def givePatientIDfromSampleName(sampleName):
	# get the information of the patient from the Name of the sample
	if sampleName.startswith('Neu'):  # if the sample name starts with
		return sampleName[3:]  # give the patient ID located in
	if sampleName.startswith('poly'):  # if the sample name starts with
		return sampleName[:5]  # give the patient ID located in
	if sampleName.startswith('Act') or sampleName.startswith('act') or sampleName.startswith(
			'Rem') or sampleName.startswith('rem'):  # if the sample name starts with
		return 'sJIA' + sampleName[3:4]  # give the patient ID located in
	if sampleName.startswith('HC'):  # if the sample name starts with
		return sampleName[:3] # give the patient ID located in
	raise Exception('no patientID for', sampleName)  # if there is no sample name, raise an error

makeMonocytesOnlyNGSbasedJunctionsFromSamFile(dataDirBase=dataDirBase, dataDirBaseBulk=dataDirBaseBulk, folderName=folderName, folderNameSamples=folderNameSamples, folderNameBatch=folderNameBatch, outputFN=outputFN)

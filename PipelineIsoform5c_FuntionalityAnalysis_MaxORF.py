#!/usr/bin/env python

#######################################
# FLAIR longest-ORF isoform assesment #
#######################################

def makeORFfastaFromACellTypeOnlyFLAIRdatabase(FileNamePredictedProtein=None, FileNamePredictedIsoforms=None, minimumORFlength=25):

	# If the FileName is not given, print an error
	if not FileNamePredictedProtein or not FileNamePredictedIsoforms:  # until there are no more lines
		raise Exception('File name not given')  # print an error in txNameOnly = TRUE

	ORFMaxORFFunctional = {}  # create a dictionary to save the functional annotation

	geneiso = {}  # create a dictionary to save the information of the genes with the isoforms
	# Open isoform annotation from FLAIR (filtered)
	INfile = open(FileNamePredictedIsoforms, 'r')
	#  open new isoform annotation
	while 1:
		aline = INfile.readline()  # Read a new line of the file
		if not aline:  # If there is not a new line (you read all the lines)
			INfile.close()  # close the connexion to the file
			break  # break the loop, go out of the while
		splitline = aline.split("\t")  # split by tab
		if splitline[2] == "transcript":
			genename=splitline[8].split(" ")[3].replace('"', "").replace('\n', "") # save the gene name
			isoname=splitline[8].split(" ")[1].replace('"', "").replace(';', "") # save the isoform name
			if genename not in geneiso:  # if the gene is not define
				geneiso[genename] = list()  # create a key with the gene and add as a value a list
			geneiso[genename].append(isoname)  # add to the list of the gene all isoform names

	isolen = {}  # create a dictionary to save the information of the length of the ORF of each isoform
	# Open isoform annotation from FLAIR (filtered)
	INfile = open(FileNamePredictedProtein, 'r')
	cnt,cntheader=0,0  # counter to 0
	#  open new isoform annotation
	while 1:
		aline = INfile.readline()  # Read a new line of the file
		if not aline:  # If there is not a new line (you read all the lines)
			INfile.close()  # close the connexion to the file
			break  # break the loop, go out of the while
		if ">" in aline:  # read header
			isoname=aline.split("\t")[0].replace(">","")  # get the isoform
			if cntheader != 0 and cntheader != cnt-2:  # the header is every two lines
				raise Exception("Not one-line fasta")  # raise an exception
			cntheader = cnt  # save the last header  # save the last header
		else:  # if not
			isolen[isoname] = len(aline.replace("\n",""))  # save the length of the ORF
		cnt = cnt + 1  # plus one

	for genename in geneiso.keys():  # enter the information of each gene
		allIsoforms=geneiso[genename]  # enter all the isoforms of a gene
		allIsoformswithORF=list(set(allIsoforms) & set(isolen.keys()))  # get only isoforms with an ORF
		if len(allIsoformswithORF) != 0:
			alllengthgene=[isolen[x] for x in allIsoformswithORF]
			MaxLengthORF=max(alllengthgene)
		for isoname in allIsoforms:  # enter all isoforms
			LenORF=0
			Value="no-predicted ORF"
			if isoname in isolen.keys():  # if it has a predicted an aa sequence
				LenORF=isolen[isoname]
				if LenORF == MaxLengthORF and MaxLengthORF > minimumORFlength:
					Value="functional"
				else:
					Value = "non-functional"
			ORFMaxORFFunctional[genename,isoname]=Value  # create the key in the dictionary

	# return the table of the functionality: only functional max ORF
	OUTfile = open(FileNamePredictedIsoforms.replace(".final.gtf",".functMaxORF.txt"), 'w')
	#  for loop to print all the information
	for genename,isoname in ORFMaxORFFunctional.keys():  # enter each key in the dictionary
		OUTfile.write(genename+'\t'+isoname+'\t'+ORFMaxORFFunctional[genename,isoname]+'\n') # write all info in a file, in an specific way
	OUTfile.close()  # close the file

makeORFfastaFromACellTypeOnlyFLAIRdatabase(FileNamePredictedProtein="FinalIsoformAnnotation/Monocyte_Isoforms_predictedprot.final.faa", FileNamePredictedIsoforms="FinalIsoformAnnotation/Monocyte_Isoforms.final.gtf")

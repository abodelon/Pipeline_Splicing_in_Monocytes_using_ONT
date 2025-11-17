#!/usr/bin/bash

echo "Process started the $(date)." # print an ok message

# MaxQuantCmd.dll --changeFolder="<new mqpar.xml>" "<new folder with fasta files>" "<new folder with raw files>" "<new folder with libraries> (only for DIA runs)" mqpar.xml
dotnet MaxQuant/MaxQuant_v2.5.1.0/bin/MaxQuantCmd.dll MaxQuant/mqpar_all_18102024.xml --changeFolder MaxQuant/mqpar_all_18102024_folder.xml MaxQuant/Reference/ MaxQuant/MassSpectrometryData/
dotnet MaxQuant/MaxQuant_v2.5.1.0/bin/MaxQuantCmd.dll MaxQuant/mqpar_all_18102024_folder.xml

echo "Process finished the $(date)." # print an ok message

sleep 900

echo "Process finished the $(date)." # print an ok message
#======================================================
# TUTORIAL
#======================================================
# git clone https://github.com/lauramilena3/Tools
# conda env create -n ViralTools -f ViralTools.yaml
# source activate ViralTools
# snakemake -j $nCores --use-conda --config input_dir=$fastaDir results_dir=$reusultDir

import os
import re
#======================================================
# Config files
#======================================================
configfile: "config.yaml"

#======================================================
# Global variables
#======================================================

INPUT_DIR =config['input_dir'].rstrip("/") 
RESULTS_DIR=config['results_dir'].rstrip("/")
dir_list = ["VIRAL_DIR"]
dir_names = [RESULTS_DIR + "/VIRAL_ID"]
dirs_dict = dict(zip(dir_list, dir_names))
SAMPLES,=glob_wildcards(INPUT_DIR + "/{sample}.fasta")
print(INPUT_DIR)
print(SAMPLES)


#======================================================
# Rules
#======================================================
 

rule all:
	input:
		expand(dirs_dict["VIRAL_DIR"]+ "/{sample}_viral_table.csv", sample=SAMPLES)

rule downloadViralTools:
	output:
		virSorter_dir=directory(config['virSorter_dir']),
		virFinder_dir=directory(config['virFinder_dir']),
	message:
		"Downloading required VirSorter and VirFinder"
	threads: 1
	shell:
		"""
		#VIRSORTER
		VS_dir="{config[virSorter_dir]}"
		echo $VS_dir
		if [ ! -d $VS_dir ]
		then
			mkdir -p tools
			cd tools
			git clone https://github.com/simroux/VirSorter.git 
			cd VirSorter/Scripts 
			make clean
			make
			cd ../../../
		fi
		#VIRFNDER
		VF_dir="{config[virFinder_dir]}"
		echo $VF_dir
   		if [ ! -d $VF_dir ]
		then
			if [ ! {config[operating_system]} == "linux" ] 
			then
				curl -OL https://raw.github.com/jessieren/VirFinder/blob/master/mac/VirFinder_1.1.tar.gz?raw=true
			else
				curl -OL https://github.com/jessieren/VirFinder/blob/master/linux/VirFinder_1.1.tar.gz?raw=true
			fi
			mkdir -p {output.virFinder_dir}
			mv VirFinder*tar.gz* {output.virFinder_dir}/VirFinder_1.1.tar.gz
		fi
		"""

rule downloadViralDB:
	output:
		virSorter_db=directory(config['virSorter_db']),
	message:
		"Downloading VirSorter database"
	threads: 1
	params:
		virSorter_db="db/VirSorter"
	shell:
		"""
		VS_db="{config[virSorter_db]}"
		echo $VS_db
		if [ ! -d $VS_db ]
		then
			curl -OL https://zenodo.org/record/1168727/files/virsorter-data-v2.tar.gz
			mkdir -p {params.virSorter_db}
			tar -xvzf virsorter-data-v2.tar.gz -C {params.virSorter_db}
			rm virsorter-data-v2.tar.gz
		fi
		"""

rule virSorter:
	input:
		representatives=INPUT_DIR + "/{sample}.fasta",
		virSorter_dir=config['virSorter_dir'],
		virSorter_db=config['virSorter_db']
	output:
		results=dirs_dict["VIRAL_DIR"] + "/{sample}_virSorter/VIRSorter_global-phage-signal.csv"
	params:
		out_folder=dirs_dict["VIRAL_DIR"] + "/{sample}_virSorter"
	message:
		"Classifing contigs with VirSorter"
	conda:
		"/viral.yaml"
	threads: 4
	shell:
		"""
		{config[virSorter_dir]}/wrapper_phage_contigs_sorter_iPlant.pl -f {input.representatives} \
			--db 2 \
			--wdir {params.out_folder} \
			--ncpu {threads} \
			--data-dir {input.virSorter_db} \
			--virome  
		"""

rule virFinder:
	input:
		representatives=INPUT_DIR + "/{sample}.fasta",
		virFinder_dir=config['virFinder_dir']
	output:
		pvalues=dirs_dict["VIRAL_DIR"] + "/{sample}_virFinder_pvalues.txt",
	params:
		virFinder_script="scripts/virfinder_wrapper.R"
	message: 
		"Scoring virus VirFinder"
	conda:
		"/viral.yaml"
	threads: 1
	shell:
		"""
		Rscript {params.virFinder_script} {input.representatives} {output.pvalues}
		"""

rule parseViralTable:
	input:
		pvalues = dirs_dict["VIRAL_DIR"] + "/{sample}_virFinder_pvalues.txt",
		categories=dirs_dict["VIRAL_DIR"] + "/{sample}_virSorter/VIRSorter_global-phage-signal.csv"
	output:
		table=dirs_dict["VIRAL_DIR"]+ "/{sample}_viral_table.csv"
	params:
		virFinder_script="scripts/virfinder_wrapper.R'",
		virFinder_dir=config['virFinder_dir']
	message: 
		"Parsing VirSorter and VirFinder results"
	threads: 1
	run:
		import pandas as pd
		results = pd.DataFrame(columns=('lenght', 'circular','type', 'VS_cat', 'VF_score', 'VF_pval'))
		#VirSorter
		VS = input.categories
		with open(VS) as fp:  
			line = fp.readline()
			cnt = 1
			while line:
				if line.startswith("#"):
					if (line.strip().split()[1].isdigit()):
						contig_type=(line.strip().split("-")[1])
				else:
					circular="N"
					contigName=line.split(",")[0].split("VIRSorter_")[1].replace(".", "_")
					category=line.split(",")[4]
					if "-circular" in contigName:
						contigName=contigName.split("-circular")[0]
						circular="Y"
					if "suggestCircular=yes" in contigName:
						circular="Y"
					results.loc[contigName, 'VS_cat'] = int(category)
					results.loc[contigName, 'circular'] = circular
					results.loc[contigName, 'type'] = contig_type
				line = fp.readline()
				cnt += 1
				
		#VirFinder
		VF = input.pvalues
		with open(VF) as fp:  
			line = fp.readline()
			cnt = 1
			while line:
				if cnt != 1:
					contigName=line.split("\t")[0].strip().replace(".", "_")
					contigLenght=line.split("\t")[1]
					contigScore=line.split("\t")[2]
					contigPval=line.split("\t")[3].split("\n")[0]
					results.loc[contigName, 'lenght'] = float(contigLenght)
					results.loc[contigName, 'VF_score'] = float(contigScore)
					results.loc[contigName, 'VF_pval'] = float(contigPval)
						#check if circular also
				line = fp.readline()
				cnt += 1	


		results.to_csv(output.table)










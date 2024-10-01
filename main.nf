#!/usr/bin/env nextflow
nextflow.enable.dsl=2



// make csv file with headers from the given input

process make_csv {
	publishDir "${params.out_dir}"
	input:
	path(fastq_input)
	output:
	path("samplelist.csv")
	
	script:
	"""
	makecsv.sh ${fastq_input}

	"""

}

//merge fastq files for each SampleName and create a merged file for each SampleNames
process merge_fastq {
	publishDir "${params.out_dir}/merged"
	label "low"
	input:
	tuple val(SampleName),path(SamplePath)
	output:
	tuple val(SampleName),path("${SampleName}.{fastq,fastq.gz}"),emit:reads

	shell:
	"""
	count=\$(ls -1 ${SamplePath}/*.gz 2>/dev/null | wc -l)
	
	
		if [[ "\${count}" != "0" ]]
		then
			cat ${SamplePath}/*.fastq.gz > ${SampleName}.fastq.gz
					
		else
			count=\$(ls -1 ${SamplePath}/*.fastq 2>/dev/null | wc -l)
			if [[ "\${count}" != "0" ]]
			then
				cat ${SamplePath}/*.fastq > ${SampleName}.fastq
				
			fi
		fi
	"""
}

//trim barcodes and adapter using porechop

process porechop {
	label "high"
	publishDir "${params.out_dir}/trimmed"
	input:
	tuple val(SampleName),path(SamplePath)
	output:
	tuple val(SampleName),path ("${SampleName}_trimmed.fastq")
	script:
	"""
	porechop -i ${SamplePath} -o ${SampleName}_trimmed.fastq
	"""
}

process dragonflye {
    label "high"
    publishDir "${params.out_dir}/Assembly",mode:"copy"
    input:
    tuple val(SampleName),path(SamplePath)
	val(medaka_model)
    output:
    val(SampleName),emit:sample
	tuple val(SampleName),path("${SampleName}_flye.fasta"),emit:assembly
	path("${SampleName}_flye-info.txt"),emit:flyeinfo
    script:
    """
    dragonflye --reads ${SamplePath} --outdir ${SampleName}_assembly --model ${medaka_model} --gsize 2.4M --nanohq --medaka 1
    # rename fasta file with samplename
    mv "${SampleName}_assembly"/flye.fasta "${SampleName}"_flye.fasta
    # rename fasta header with samplename
    sed -i 's/contig/${SampleName}_contig/g' "${SampleName}_flye.fasta"
     # rename flyeinfo file and contnents
    mv "${SampleName}_assembly"/flye-info.txt "${SampleName}"_flye-info.txt
    sed -i 's/contig/${SampleName}_contig/g' "${SampleName}_flye-info.txt"
    """
}


process busco {
    label "low"
    publishDir "${params.out_dir}/busco",mode:"copy"
    input:
    tuple val(SampleName),path(cons)
    output:
    path ("${SampleName}_busco.txt")
    script:

    """
    busco -i ${cons} -m genome -l bacteria_odb10 -o ${SampleName}_busco_results
	mv ${SampleName}_busco_results/*.txt ${SampleName}_busco.txt

    """
}

process mlst {
	publishDir "${params.out_dir}/mlst/",mode:"copy"
	label "low"
	input:
	tuple val(SampleName),path(consensus)
	output:
	path("${SampleName}_MLST.csv")
	script:
	"""
	mlst ${consensus} > ${SampleName}_MLST.csv
	sed -i 's,_flye.fasta,,g' ${SampleName}_MLST.csv
	"""
}

process abricate{
	publishDir "${params.out_dir}/abricate/",mode:"copy"
	label "low"
	input:
	tuple val(SampleName),path(consensus)
	output:
	path("${SampleName}_sero.csv"),emit:sero
	path("${SampleName}_vf.csv"),emit:vif
	path("${SampleName}_AMR.csv"),emit:AMR
	
	script:
	"""
	abricate --db ecoh ${consensus} 1> ${SampleName}_sero.csv
	sed -i 's,_flye.fasta,,g' ${SampleName}_sero.csv
	abricate --db ecoli_vf  ${consensus} 1> ${SampleName}_vf.csv
	sed -i 's,_flye.fasta,,g' ${SampleName}_vf.csv
	abricate --db card ${consensus} 1> ${SampleName}_AMR.csv
	sed -i 's,_flye.fasta,,g' ${SampleName}_AMR.csv
	
	"""
	
}

   

process make_limsfile {
	label "low"
	publishDir "${params.out_dir}/LIMS",mode:"copy"
	input:
	path (serotyping_results)
	path (vf_results)
	path (amr_results)
	path (mlst_results)
	path (software_version)
	output:
	path("*_LIMS_file.csv")
	path("sero_file.csv"),emit:sero
	path("MLST_file.csv"),emit:mlst
	
	script:
	"""
	LIMS_file.sh
	
	
	awk 'FNR==1 && NR!=1 { while (/^#F/) getline; } 1 {print}' ${mlst_results} > MLST_file.csv
	# add header to mlst file
	sed -i \$'1 i\\\nSAMPLE\tSCHEME\tST\tadk\tfumC\tgyrB\ticd\tmdh\tpurA\trecA' MLST_file.csv

	"""
}

process make_report {
	label "low"
	publishDir "${params.out_dir}/",mode:"copy"
	input:
	path(rmdfile)
	path(lims)
	path(serofile)
	path (mlstfile)
	path(busco)
	path (samplelist)
	path (vffiles)
	path (amrfiles)
	output:
	path("Ecoli_report.html")

	script:

	"""
	
	cp ${rmdfile} rmdfile_copy.Rmd
	cp ${samplelist} samples.csv
	cp ${serofile} serofile.csv
	cp ${mlstfile} mlstfile.csv

	Rscript -e 'rmarkdown::render(input="rmdfile_copy.Rmd",params=list(sero="serofile.csv",mlst="mlstfile.csv",csv="samples.csv"),output_file="Ecoli_report.html")'
	
	
	"""

}





workflow {
    data=Channel
	.fromPath(params.input)
	merge_fastq(make_csv(data).splitCsv(header:true).map { row-> tuple(row.SampleName,row.SamplePath)})
	// Merge fastq files for each sample

	// based on the optional argument trim barcodes using porechop and assemble using dragonflye
    if (params.trim_barcodes){
		porechop(merge_fastq.out)
		dragonflye(porechop.out,params.medaka_model) 
	} else {
        dragonflye(merge_fastq.out,params.medaka_model)           
    }
	versionfile=file("${baseDir}/software_version.csv")
	//checking completeness of assembly
    busco(dragonflye.out.assembly)
	//mlst
    mlst (dragonflye.out.assembly)
	
    abricate (dragonflye.out.assembly)
	
	
    //make lims file
    make_limsfile (abricate.out.sero.collect(),abricate.out.vif.collect(),abricate.out.AMR.collect(),mlst.out.collect(),versionfile)
	//report generation

	rmd_file=file("${baseDir}/Ecoli_report.Rmd")
	make_report (rmd_file,make_limsfile.out,busco.out.collect(),make_csv.out,abricate.out.vif.collect(),abricate.out.AMR.collect())


}


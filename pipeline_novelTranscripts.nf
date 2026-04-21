#!/usr/bin/env nextflow

/*
 * Pipeline parameters
 */


//Output directories


/*
 * Processes
 */


//suppa generate exon skipping events from model gtf
process suppa_ge {

        publishDir params.outdir

        container params.spliceCont

	input:
            path models_gtf

	output:
	    tuple path("${models_gtf}"), path("transcript_model_as_events_SE_strict.ioe")

	script:
	"""
        python /usr/local/suppa/suppa.py generateEvents -i '$models_gtf' \
                                -o transcript_model_as_events \
                                -f ioe \
                                -e SE
	"""
}

//suppa calculate psi (relative abundance) for each transcript
process suppa_isoPsi {

        publishDir params.outdir

        container params.spliceCont

        input:
            tuple path(models_gtf), path(se_ioe)
	    path model_tpms

        output:
            tuple path("${models_gtf}"), path("${se_ioe}"), path("transcript_model_isoform.psi")

        script:
        """
        python /usr/local/suppa/suppa.py psiPerIsoform -g '$models_gtf' \
                               -e '$model_tpms' \
                               -o transcript_model
        """
}

//convert se ioe file to bed file of cassette exons
process ioe2bed {

        publishDir params.outdir

        container params.rCont

        input:
            tuple path(models_gtf), path(se_ioe), path(iso_psi)

        output:
            tuple path("${models_gtf}"), path("cassetteExon.bed"), path("${iso_psi}")

        script:
        """
        Rscript '$projectDir'/scripts/ioe2bed.R -e '$se_ioe'
        """
}

//overlap cassette exon bed with gencode ref gtf exons using bedtools intersect -wao (null overlaps reported) 
process wao_intersect {

        publishDir params.outdir

        container params.spliceCont

        input:
            tuple path(models_gtf), path(ce_bed), path(iso_psi)
            path refgtf

        output:
            tuple path("${models_gtf}"), path("ce_vs_ref.txt"), path("${iso_psi}")

        script:
        """
        bedtools intersect -wao -s -a '$ce_bed' -b <(zcat '$refgtf' | awk '\$3=="exon"') > ce_vs_ref.txt
        """
}

//putative cassette exon QC, associate with transcripts and write final output tables
process get_novel_transcripts {

        publishDir params.outdir

        container params.rCont

        input:
            tuple path(models_gtf), path(ce_vs_ref), path(iso_psi)
            path model_tpms

        output:
            tuple path("novelTranscripts.tsv"), path("novelCassetteExons.tsv")

        script:
        """
        Rscript '$projectDir'/scripts/getNovelTranscripts.R -c '$ce_vs_ref' \
                                                            -m '$models_gtf' \
                                                            -p '$iso_psi' \
                                                            -t '$model_tpms'
        """
}

/*
 * Workflow
 */


workflow {

    //create input channel
    models_ch = channel.fromPath(params.models)
    
    suppa_ge(models_ch)

    suppa_isoPsi(suppa_ge.out, file(params.model_tpms))

    ioe2bed(suppa_isoPsi.out)

    wao_intersect(ioe2bed.out), file(params.refgtf)

    get_novel_transcripts(wao_intersect.out, file(params.model_tpms))
}


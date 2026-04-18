version 1.0

workflow MosaicPostProcessing {
    input {
        File proband_hc_vcf
        File proband_hc_idx
        File proband_mutect_vcf
        File proband_mutect_idx
        File mother_hc_vcf
        File mother_hc_idx
        File father_hc_vcf
        File father_hc_idx
        
        File proband_bam
        File proband_bai
        File mother_bam
        File mother_bai
        File father_bam
        File father_bai
        
        File ref_fasta
        File ref_index
        
        File? centromere_bed
        File? repetitive_bed
        Float min_tlod = 3.0
        
        File postprocess_script
    }

    # Path A: HC Filtration
    call MosaicPostProcess as Filter_HC {
        input:
            proband_vcf = proband_hc_vcf,
            proband_idx = proband_hc_idx,
            mother_vcf = mother_hc_vcf,
            mother_idx = mother_hc_idx,
            father_vcf = father_hc_vcf,
            father_idx = father_hc_idx,
            proband_bam = proband_bam,
            proband_bai = proband_bai,
            mother_bam = mother_bam,
            mother_bai = mother_bai,
            father_bam = father_bam,
            father_bai = father_bai,
            ref_fasta = ref_fasta,
            ref_index = ref_index,
            centromere_bed = centromere_bed,
            repetitive_bed = repetitive_bed,
            min_tlod = min_tlod,
            script = postprocess_script,
            output_prefix = "mosaic_candidates_hc"
    }

    # Path B: Mutect Filtration
    call MosaicPostProcess as Filter_Mutect {
        input:
            proband_vcf = proband_mutect_vcf,
            proband_idx = proband_mutect_idx,
            mother_vcf = mother_hc_vcf,
            mother_idx = mother_hc_idx,
            father_vcf = father_hc_vcf,
            father_idx = father_hc_idx,
            proband_bam = proband_bam,
            proband_bai = proband_bai,
            mother_bam = mother_bam,
            mother_bai = mother_bai,
            father_bam = father_bam,
            father_bai = father_bai,
            ref_fasta = ref_fasta,
            ref_index = ref_index,
            centromere_bed = centromere_bed,
            repetitive_bed = repetitive_bed,
            min_tlod = min_tlod,
            script = postprocess_script,
            output_prefix = "mosaic_candidates_mutect"
    }

    output {
        File candidates_hc_tsv = Filter_HC.candidates_tsv
        File candidates_mutect_tsv = Filter_Mutect.candidates_tsv
        File report_hc = Filter_HC.report
        File report_mutect = Filter_Mutect.report
    }
}

task MosaicPostProcess {
    input {
        File proband_vcf
        File proband_idx
        File mother_vcf
        File mother_idx
        File father_vcf
        File father_idx
        
        File proband_bam
        File proband_bai
        File mother_bam
        File mother_bai
        File father_bam
        File father_bai
        
        File ref_fasta
        File ref_index
        
        File? centromere_bed
        File? repetitive_bed
        Float min_tlod = 3.0
        
        File script
        String output_prefix
    }
    
    command <<<
        set -e
        # Dependency check
        if ! command -v bcftools &> /dev/null; then
            apt-get update && apt-get install -y bcftools
        fi
        if ! command -v pip &> /dev/null; then
             apt-get update && apt-get install -y python3-pip
        fi
        pip install pysam
        
        python3 ~{script} \
            --proband-vcf ~{proband_vcf} \
            --mother-vcf ~{mother_vcf} \
            --father-vcf ~{father_vcf} \
            --proband-bam ~{proband_bam} \
            --mother-bam ~{mother_bam} \
            --father-bam ~{father_bam} \
            --ref-fasta ~{ref_fasta} \
            ~{"--centromere-bed " + centromere_bed} \
            ~{"--repetitive-bed " + repetitive_bed} \
            --min-tlod ~{min_tlod} \
            --output-prefix ~{output_prefix} \
            --work-dir "work"
    >>>

    runtime {
        docker: "python:3.9-bullseye"
        memory: "4 GB"
        cpu: 1
        disks: "local-disk 20 HDD"
    }
    
    output {
        File candidates_tsv = "~{output_prefix}.mosaic_candidates.tsv"
        File report = "~{output_prefix}_report.txt"
    }
}

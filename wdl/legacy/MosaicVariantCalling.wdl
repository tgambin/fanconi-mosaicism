version 1.0

workflow MosaicVariantCalling {
    input {
        File proband_bam
        File proband_bai
        File mother_bam
        File mother_bai
        File father_bam
        File father_bai
        
        File ref_fasta
        File ref_index
        File ref_dict
        
        String proband_name = "proband"
        String mother_name = "mother"
        String father_name = "father"
        
        File? intervals
        
        # Mutect2 optional resources
        File? pon
        File? pon_index
        File? germline_resource
        File? germline_resource_index
        
        Int cpu = 2
        String memory = "8 GB"
    }

    call HaplotypeCallerDefault as HC_Proband {
        input:
            input_bam = proband_bam,
            input_bai = proband_bai,
            ref_fasta = ref_fasta,
            ref_index = ref_index,
            ref_dict = ref_dict,
            sample_name = proband_name,
            intervals = intervals,
            cpu = cpu,
            memory = memory
    }

    call HaplotypeCallerDefault as HC_Mother {
        input:
            input_bam = mother_bam,
            input_bai = mother_bai,
            ref_fasta = ref_fasta,
            ref_index = ref_index,
            ref_dict = ref_dict,
            sample_name = mother_name,
            intervals = intervals,
            cpu = cpu,
            memory = memory
    }

    call HaplotypeCallerDefault as HC_Father {
        input:
            input_bam = father_bam,
            input_bai = father_bai,
            ref_fasta = ref_fasta,
            ref_index = ref_index,
            ref_dict = ref_dict,
            sample_name = father_name,
            intervals = intervals,
            cpu = cpu,
            memory = memory
    }

    call Mutect2TumorOnly {
        input:
            input_bam = proband_bam,
            input_bai = proband_bai,
            ref_fasta = ref_fasta,
            ref_index = ref_index,
            ref_dict = ref_dict,
            sample_name = proband_name,
            intervals = intervals,
            pon = pon,
            pon_index = pon_index,
            germline_resource = germline_resource,
            germline_resource_index = germline_resource_index,
            cpu = cpu,
            memory = memory
    }

    output {
        File hc_proband_vcf = HC_Proband.output_vcf
        File hc_proband_idx = HC_Proband.output_vcf_index
        File hc_mother_vcf = HC_Mother.output_vcf
        File hc_mother_idx = HC_Mother.output_vcf_index
        File hc_father_vcf = HC_Father.output_vcf
        File hc_father_idx = HC_Father.output_vcf_index
        File mutect_proband_vcf = Mutect2TumorOnly.output_vcf
        File mutect_proband_idx = Mutect2TumorOnly.output_vcf_index
    }
}

task HaplotypeCallerDefault {
    input {
        File input_bam
        File input_bai
        File ref_fasta
        File ref_index
        File ref_dict
        String sample_name
        File? intervals
        Int ploidy = 2
        
        Int cpu
        String memory
    }

    command {
        gatk HaplotypeCaller \
            -R ~{ref_fasta} \
            -I ~{input_bam} \
            -O ~{sample_name}.hc.ploidy~{ploidy}.vcf.gz \
            -ploidy ~{ploidy} \
            --native-pair-hmm-threads ~{cpu} \
            ~{"-L " + intervals}
    }

    runtime {
        docker: "broadinstitute/gatk:latest"
        memory: memory
        cpu: cpu
        disks: "local-disk 50 HDD"
    }

    output {
        File output_vcf = "~{sample_name}.hc.ploidy~{ploidy}.vcf.gz"
        File output_vcf_index = "~{sample_name}.hc.ploidy~{ploidy}.vcf.gz.tbi"
    }
}

task Mutect2TumorOnly {
    input {
        File input_bam
        File input_bai
        File ref_fasta
        File ref_index
        File ref_dict
        String sample_name
        File? intervals
        
        File? pon
        File? pon_index
        File? germline_resource
        File? germline_resource_index
        
        Int cpu
        String memory
    }

    command {
        gatk Mutect2 \
            -R ~{ref_fasta} \
            -I ~{input_bam} \
            -O ~{sample_name}.mutect2.vcf.gz \
            --native-pair-hmm-threads ~{cpu} \
            ~{"-L " + intervals} \
            ~{"--panel-of-normals " + pon} \
            ~{"--germline-resource " + germline_resource}
    }

    runtime {
        docker: "broadinstitute/gatk:latest"
        memory: memory
        cpu: cpu
        disks: "local-disk 50 HDD"
    }

    output {
        File output_vcf = "~{sample_name}.mutect2.vcf.gz"
        File output_vcf_index = "~{sample_name}.mutect2.vcf.gz.tbi"
    }
}

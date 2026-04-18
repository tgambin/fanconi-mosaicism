version 1.0

workflow Mutect2Standalone {
    input {
        File proband_bam
        File proband_bai
        File ref_fasta
        File ref_index
        File ref_dict
        
        String sample_name = "proband"
        File? intervals
        
        # Mutect2 optional resources
        File? pon
        File? pon_index
        File? germline_resource
        File? germline_resource_index
        
        Int cpu = 2
        String memory = "8 GB"
        Int disk_size = 50
    }

    call Mutect2 {
        input:
            input_bam = proband_bam,
            input_bai = proband_bai,
            ref_fasta = ref_fasta,
            ref_index = ref_index,
            ref_dict = ref_dict,
            sample_name = sample_name,
            intervals = intervals,
            pon = pon,
            pon_index = pon_index,
            germline_resource = germline_resource,
            germline_resource_index = germline_resource_index,
            cpu = cpu,
            memory = memory,
            disk_size = disk_size
    }

    call FilterMutectCalls {
        input:
            unfiltered_vcf = Mutect2.output_vcf,
            unfiltered_vcf_idx = Mutect2.output_vcf_index,
            stats = Mutect2.output_stats,
            ref_fasta = ref_fasta,
            ref_index = ref_index,
            ref_dict = ref_dict,
            sample_name = sample_name,
            memory = memory,
            disk_size = disk_size
    }

    output {
        File raw_vcf = Mutect2.output_vcf
        File raw_vcf_idx = Mutect2.output_vcf_index
        File stats = Mutect2.output_stats
        File filtered_vcf = FilterMutectCalls.filtered_vcf
        File filtered_vcf_idx = FilterMutectCalls.filtered_vcf_index
    }
}

task Mutect2 {
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
        Int disk_size
    }

    command {
        gatk Mutect2 \
            -R ~{ref_fasta} \
            -I ~{input_bam} \
            -O ~{sample_name}.raw.mutect2.vcf.gz \
            --native-pair-hmm-threads ~{cpu} \
            ~{"-L " + intervals} \
            ~{"--panel-of-normals " + pon} \
            ~{"--germline-resource " + germline_resource}
    }

    runtime {
        docker: "broadinstitute/gatk:latest"
        memory: memory
        cpu: cpu
        disks: "local-disk " + disk_size + " HDD"
    }

    output {
        File output_vcf = "~{sample_name}.raw.mutect2.vcf.gz"
        File output_vcf_index = "~{sample_name}.raw.mutect2.vcf.gz.tbi"
        File output_stats = "~{sample_name}.raw.mutect2.vcf.gz.stats"
    }
}

task FilterMutectCalls {
    input {
        File unfiltered_vcf
        File unfiltered_vcf_idx
        File stats
        File ref_fasta
        File ref_index
        File ref_dict
        String sample_name
        
        String memory
        Int disk_size
    }

    command {
        gatk FilterMutectCalls \
            -R ~{ref_fasta} \
            -V ~{unfiltered_vcf} \
            --stats ~{stats} \
            -O ~{sample_name}.filtered.vcf.gz
    }

    runtime {
        docker: "broadinstitute/gatk:latest"
        memory: memory
        cpu: 1
        disks: "local-disk " + disk_size + " HDD"
    }

    output {
        File filtered_vcf = "~{sample_name}.filtered.vcf.gz"
        File filtered_vcf_index = "~{sample_name}.filtered.vcf.gz.tbi"
    }
}

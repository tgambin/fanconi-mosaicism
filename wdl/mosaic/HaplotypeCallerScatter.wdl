version 1.0

workflow HaplotypeCallerScatter {
    input {
        String proband_bam
        String proband_bai
        File ref_fasta
        File ref_index
        File ref_dict

        String sample_name = "proband"
        File intervals
        Int scatter_count = 10
        Int sample_ploidy = 50

        Int cpu_per_scatter = 2
        String memory_per_scatter = "4 GB"
        Int disk_size = 50
        Int disk_size_per_scatter = 50
        Int preemptible_attempts = 3
    }

    call SplitIntervals {
        input:
            ref_fasta = ref_fasta,
            ref_index = ref_index,
            ref_dict = ref_dict,
            intervals = intervals,
            scatter_count = scatter_count,
            preemptible_attempts = preemptible_attempts
    }

    scatter (interval_file in SplitIntervals.interval_files) {
        call HaplotypeCaller {
            input:
                input_bam = proband_bam,
                input_bai = proband_bai,
                ref_fasta = ref_fasta,
                ref_index = ref_index,
                ref_dict = ref_dict,
                sample_name = sample_name,
                intervals = interval_file,
                sample_ploidy = sample_ploidy,
                cpu = cpu_per_scatter,
                memory = memory_per_scatter,
                disk_size = disk_size_per_scatter,
                preemptible_attempts = preemptible_attempts
        }
    }

    call MergeVCFs {
        input:
            vcfs = HaplotypeCaller.output_vcf,
            vcf_indices = HaplotypeCaller.output_vcf_index,
            output_name = sample_name + ".hc_ploidy" + sample_ploidy,
            disk_size = disk_size,
            preemptible_attempts = preemptible_attempts
    }

    output {
        File output_vcf = MergeVCFs.output_vcf
        File output_vcf_index = MergeVCFs.output_vcf_index
    }
}

task SplitIntervals {
    input {
        File ref_fasta
        File ref_index
        File ref_dict
        File intervals
        Int scatter_count
        String memory = "4 GB"
        Int disk_size = 20
        Int preemptible_attempts = 3
    }

    command <<<
        mkdir intervals_dir
        gatk SplitIntervals \
            -R ~{ref_fasta} \
            -L ~{intervals} \
            --scatter-count ~{scatter_count} \
            --subdivision-mode INTERVAL_SUBDIVISION \
            -O intervals_dir
    >>>

    output {
        Array[File] interval_files = glob("intervals_dir/*.interval_list")
    }

    runtime {
        docker: "broadinstitute/gatk:latest"
        memory: memory
        cpu: 1
        disks: "local-disk " + disk_size + " HDD"
        preemptible: preemptible_attempts
        maxRetries: 3
        bootDiskSizeGb: 20
    }
}

task HaplotypeCaller {
    input {
        String input_bam
        String input_bai
        File ref_fasta
        File ref_index
        File ref_dict
        String sample_name
        File intervals
        Int sample_ploidy

        Int cpu
        String memory
        Int disk_size
        Int preemptible_attempts = 3
    }

    command {
        gatk HaplotypeCaller \
            -R ~{ref_fasta} \
            -I ~{input_bam} \
            --read-index ~{input_bai} \
            -O ~{sample_name}.hc.vcf.gz \
            --sample-ploidy ~{sample_ploidy} \
            --native-pair-hmm-threads ~{cpu} \
            -L ~{intervals}
    }

    runtime {
        docker: "broadinstitute/gatk:latest"
        memory: memory
        cpu: cpu
        disks: "local-disk " + disk_size + " HDD"
        preemptible: preemptible_attempts
        maxRetries: 3
        bootDiskSizeGb: 20
    }

    output {
        File output_vcf = "~{sample_name}.hc.vcf.gz"
        File output_vcf_index = "~{sample_name}.hc.vcf.gz.tbi"
    }
}

task MergeVCFs {
    input {
        Array[File] vcfs
        Array[File] vcf_indices
        String output_name
        String memory = "8 GB"
        Int disk_size = 100
        Int preemptible_attempts = 3
    }

    command <<<
        gatk MergeVcfs \
            -I ~{sep=" -I " vcfs} \
            -O ~{output_name}.vcf.gz
    >>>

    output {
        File output_vcf = "~{output_name}.vcf.gz"
        File output_vcf_index = "~{output_name}.vcf.gz.tbi"
    }

    runtime {
        docker: "broadinstitute/gatk:latest"
        memory: memory
        cpu: 1
        disks: "local-disk " + disk_size + " HDD"
        preemptible: preemptible_attempts
        maxRetries: 3
        bootDiskSizeGb: 20
    }
}

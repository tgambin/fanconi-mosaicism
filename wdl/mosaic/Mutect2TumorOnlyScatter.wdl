version 1.0

workflow Mutect2TumorOnlyScatter {
    input {
        String proband_bam
        String proband_bai
        File ref_fasta
        File ref_index
        File ref_dict
        
        String sample_name = "proband"
        File intervals
        Int scatter_count = 10
        
        # Mutect2 optional resources
        File? pon
        File? pon_index
        File? germline_resource
        File? germline_resource_index
        
        Int cpu_per_scatter = 2
        String memory_per_scatter = "4 GB"
        Int disk_size = 50
        Int disk_size_per_scatter = 50
    }

    call SplitIntervals {
        input:
            ref_fasta = ref_fasta,
            ref_index = ref_index,
            ref_dict = ref_dict,
            intervals = intervals,
            scatter_count = scatter_count
    }

    scatter (interval_file in SplitIntervals.interval_files) {
        call Mutect2 {
            input:
                input_bam = proband_bam,
                input_bai = proband_bai,
                ref_fasta = ref_fasta,
                ref_index = ref_index,
                ref_dict = ref_dict,
                sample_name = sample_name,
                intervals = interval_file,
                pon = pon,
                pon_index = pon_index,
                germline_resource = germline_resource,
                germline_resource_index = germline_resource_index,
                cpu = cpu_per_scatter,
                memory = memory_per_scatter,
                disk_size = disk_size_per_scatter
        }
    }

    call MergeVCFs {
        input:
            vcfs = Mutect2.output_vcf,
            vcf_indices = Mutect2.output_vcf_index,
            output_name = sample_name + ".raw.mutect2"
    }

    call MergeMutectStats {
        input:
            stats = Mutect2.output_stats,
            output_name = sample_name + ".raw.mutect2"
    }

    call FilterMutectCalls {
        input:
            unfiltered_vcf = MergeVCFs.output_vcf,
            unfiltered_vcf_idx = MergeVCFs.output_vcf_index,
            stats = MergeMutectStats.merged_stats,
            ref_fasta = ref_fasta,
            ref_index = ref_index,
            ref_dict = ref_dict,
            sample_name = sample_name,
            memory = "8 GB",
            disk_size = disk_size
    }

    output {
        File raw_vcf = MergeVCFs.output_vcf
        File raw_vcf_idx = MergeVCFs.output_vcf_index
        File stats = MergeMutectStats.merged_stats
        File filtered_vcf = FilterMutectCalls.filtered_vcf
        File filtered_vcf_idx = FilterMutectCalls.filtered_vcf_index
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
        preemptible: 5
        maxRetries: 5
        bootDiskSizeGb: 20
    }
}

task Mutect2 {
    input {
        String input_bam
        String input_bai
        File ref_fasta
        File ref_index
        File ref_dict
        String sample_name
        File intervals
        
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
            --read-index ~{input_bai} \
            -O ~{sample_name}.raw.mutect2.vcf.gz \
            --native-pair-hmm-threads ~{cpu} \
            -L ~{intervals} \
            ~{"--panel-of-normals " + pon} \
            ~{"--germline-resource " + germline_resource}
    }

    runtime {
        docker: "broadinstitute/gatk:latest"
        memory: memory
        cpu: cpu
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 5
        maxRetries: 5
        bootDiskSizeGb: 20
    }

    output {
        File output_vcf = "~{sample_name}.raw.mutect2.vcf.gz"
        File output_vcf_index = "~{sample_name}.raw.mutect2.vcf.gz.tbi"
        File output_stats = "~{sample_name}.raw.mutect2.vcf.gz.stats"
    }
}

task MergeVCFs {
    input {
        Array[File] vcfs
        Array[File] vcf_indices
        String output_name
        String memory = "8 GB"
        Int disk_size = 100
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
        preemptible: 5
        maxRetries: 5
        bootDiskSizeGb: 20
    }
}

task MergeMutectStats {
    input {
        Array[File] stats
        String output_name
        String memory = "4 GB"
        Int disk_size = 20
    }

    command <<<
        gatk MergeMutectStats \
            -stats ~{sep=" -stats " stats} \
            -O ~{output_name}.stats
    >>>

    output {
        File merged_stats = "~{output_name}.stats"
    }

    runtime {
        docker: "broadinstitute/gatk:latest"
        memory: memory
        cpu: 1
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 5
        maxRetries: 5
        bootDiskSizeGb: 20
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
        preemptible: 5
        maxRetries: 5
        bootDiskSizeGb: 20
    }

    output {
        File filtered_vcf = "~{sample_name}.filtered.vcf.gz"
        File filtered_vcf_index = "~{sample_name}.filtered.vcf.gz.tbi"
    }
}

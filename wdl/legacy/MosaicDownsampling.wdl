version 1.0

workflow MosaicDownsampling {
    input {
        File input_bam
        File input_bai
        Float probability
        String output_prefix
    }

    call Downsample {
        input:
            bam = input_bam,
            bai = input_bai,
            p = probability,
            prefix = output_prefix
    }

    output {
        File downsampled_bam = Downsample.out_bam
        File downsampled_bai = Downsample.out_bai
    }
}

task Downsample {
    input {
        File bam
        File bai
        Float p
        String prefix
    }

    command {
        gatk DownsampleSam \
            -I ~{bam} \
            -O ~{prefix}.bam \
            -P ~{p} \
            --CREATE_INDEX true \
            --VALIDATION_STRINGENCY SILENT
    }

    output {
        File out_bam = "~{prefix}.bam"
        File out_bai = "~{prefix}.bai"
    }

    runtime {
        docker: "broadinstitute/gatk:4.6.0.0"
        memory: "4 GB"
        cpu: 1
    }
}

version 1.0

## AddReadGroups: Add @RG header to BAM if missing
##
## Checks if BAM already has @RG. If yes, outputs the original BAM unchanged.
## If not, adds @RG via gatk AddOrReplaceReadGroups.
## Safe to run on any BAM — idempotent.

workflow AddReadGroups {
    input {
        File input_bam
        File input_bai
        String sample_name
        String library = "lib1"
        String platform = "ILLUMINA"
        String platform_unit = "unit1"
    }

    call FixReadGroups {
        input:
            input_bam = input_bam,
            input_bai = input_bai,
            sample_name = sample_name,
            library = library,
            platform = platform,
            platform_unit = platform_unit
    }

    output {
        File output_bam = FixReadGroups.output_bam
        File output_bai = FixReadGroups.output_bai
    }
}

task FixReadGroups {
    input {
        File input_bam
        File input_bai
        String sample_name
        String library
        String platform
        String platform_unit
    }

    Int disk_size = ceil(size(input_bam, "GB") * 2.5) + 50

    command <<<
        set -euo pipefail

        # Check if @RG already exists in BAM header
        RG_COUNT=$(samtools view -H ~{input_bam} | grep -c "^@RG" || true)

        if [ "$RG_COUNT" -gt 0 ]; then
            echo "@RG already present ($RG_COUNT lines) — passing through unchanged"
            cp ~{input_bam} ~{sample_name}.rg.bam
            cp ~{input_bai} ~{sample_name}.rg.bam.bai
        else
            echo "No @RG found — adding RG header + read tags via samtools addreplacerg"
            samtools addreplacerg -@ 4 \
                -r "ID:~{sample_name}\tSM:~{sample_name}\tLB:~{library}\tPL:~{platform}\tPU:~{platform_unit}" \
                -o ~{sample_name}.rg.bam ~{input_bam}
            samtools index -@ 4 ~{sample_name}.rg.bam
        fi
    >>>

    runtime {
        docker: "staphb/samtools:1.19"
        cpu: 4
        memory: "16 GB"
        disks: "local-disk " + disk_size + " SSD"
        preemptible: 1
    }

    output {
        File output_bam = "~{sample_name}.rg.bam"
        File output_bai = "~{sample_name}.rg.bam.bai"
    }
}

version 1.0

workflow MosdepthCoverage {
    input {
        File bam
        File bai
        String sample_name
        File? coverage_bed
    }

    call RunMosdepth {
        input:
            bam = bam,
            bai = bai,
            sample = sample_name,
            bed = coverage_bed
    }

    output {
        File summary = RunMosdepth.summary
        File global_dist = RunMosdepth.global_dist
        File? regions_bed = RunMosdepth.regions_bed
        File? regions_mean = RunMosdepth.regions_mean
    }
}

task RunMosdepth {
    input {
        File bam
        File bai
        String sample
        File? bed
    }

    command <<<
        mosdepth --fast-mode --no-per-base ~{"--by " + bed} ~{sample} ~{bam}
        
        # If BED was provided, extract the mean coverage for the first region
        if [ -n "~{bed}" ]; then
            # The regions.bed.gz contains: chrom start end mean_coverage
            # We unzip it and extract the 4th column (mean)
            zcat ~{sample}.regions.bed.gz | head -1 | cut -f4 > ~{sample}.mean_coverage.txt
        fi
    >>>

    output {
        File summary = "~{sample}.mosdepth.summary.txt"
        File global_dist = "~{sample}.mosdepth.global.dist.txt"
        File? regions_bed = "~{sample}.regions.bed.gz"
        File? regions_mean = "~{sample}.mean_coverage.txt"
    }

    runtime {
        docker: "quay.io/biocontainers/mosdepth:0.3.3--h37c5b7d_2"
        memory: "2 GB"
        cpu: 1
    }
}

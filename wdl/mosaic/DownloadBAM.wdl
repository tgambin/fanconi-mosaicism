version 1.0

## DownloadBAM: Download a BAM+BAI from HTTPS URL to GCS via Terra
##
## Simple wrapper — runs wget on a Terra VM, outputs BAM+BAI.
## Cromwell delocalizes outputs to GCS automatically.
## Designed for GIAB 300x BAMs from NCBI (~560 GB each).

workflow DownloadBAM {
    input {
        String bam_url
        String bai_url
        String sample_name
        Int disk_size = 700
    }

    call Download {
        input:
            bam_url = bam_url,
            bai_url = bai_url,
            sample_name = sample_name,
            disk_size = disk_size
    }

    output {
        File bam = Download.bam
        File bai = Download.bai
    }
}

task Download {
    input {
        String bam_url
        String bai_url
        String sample_name
        Int disk_size
    }

    command <<<
        set -euo pipefail

        apt-get update && apt-get install -y curl

        echo "Downloading ~{sample_name} BAM from ~{bam_url}"
        echo "Start: $(date)"

        curl -sSL -o ~{sample_name}.bam "~{bam_url}"
        echo "BAM done: $(du -h ~{sample_name}.bam | cut -f1) — $(date)"

        curl -sSL -o ~{sample_name}.bam.bai "~{bai_url}"
        echo "BAI done — $(date)"
    >>>

    runtime {
        docker: "ubuntu:22.04"
        memory: "2 GB"
        cpu: 2
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
        maxRetries: 1
    }

    output {
        File bam = "~{sample_name}.bam"
        File bai = "~{sample_name}.bam.bai"
    }
}

version 1.0

## IGVTrio: BAM slicing + IGV headless screenshots for a trio
## Single-file WDL — no imports needed for Terra.
##
## Steps:
##   1. GenerateRegions: parse TSV → samtools region args
##   2. SliceBam × 3: stream BAM slices from GCS (proband, mother, father)
##   3. CreateScreenshots: IGV headless batch on mini-BAMs → tar.gz

workflow IGVTrio_Final {
    input {
        File variants_tsv
        String proband_bam
        File proband_bai
        String mother_bam
        File mother_bai
        String father_bam
        File father_bai
        File ref_fasta
        File ref_index

        String output_name = "igv_screenshots"
        Int padding = 200
        Boolean is_terra = true
        String igv_docker_image = "tgambin/igv-headless:latest"
    }

    call GenerateRegions {
        input:
            tsv = variants_tsv,
            pad = padding
    }

    call SliceBam as SliceProband {
        input:
            bam = proband_bam,
            bai = proband_bai,
            regions_args = GenerateRegions.regions_string,
            sample_name = "proband",
            is_terra = is_terra
    }

    call SliceBam as SliceMother {
        input:
            bam = mother_bam,
            bai = mother_bai,
            regions_args = GenerateRegions.regions_string,
            sample_name = "mother",
            is_terra = is_terra
    }

    call SliceBam as SliceFather {
        input:
            bam = father_bam,
            bai = father_bai,
            regions_args = GenerateRegions.regions_string,
            sample_name = "father",
            is_terra = is_terra
    }

    call CreateScreenshots {
        input:
            variants_tsv = variants_tsv,
            p_bam = SliceProband.minibam,
            p_bai = SliceProband.minibai,
            m_bam = SliceMother.minibam,
            m_bai = SliceMother.minibai,
            f_bam = SliceFather.minibam,
            f_bai = SliceFather.minibai,
            ref_fasta = ref_fasta,
            ref_index = ref_index,
            out_name = output_name,
            pad = padding,
            docker_image = igv_docker_image,
            is_terra = is_terra
    }

    output {
        File screenshots_tar = CreateScreenshots.tar_gz
    }
}

# --- Tasks ---

task GenerateRegions {
    input {
        File tsv
        Int pad
    }

    command <<<
        awk -v p=~{pad} '
        BEGIN { ORS=" " }
        {
            sub(/\r$/, "", $0);
            if ($2 !~ /^[0-9]+$/) next;
            chrom = $1;
            if (index(tolower(chrom), "chr") != 1) chrom = "chr" chrom;
            start = $2 - p;
            if (start < 1) start = 1;
            end = $2 + p;
            print chrom ":" start "-" end
        }' "~{tsv}" > regions.txt
    >>>

    runtime {
        docker: "ubuntu:22.04"
        memory: "1 GB"
        cpu: 1
        disks: "local-disk 10 HDD"
        preemptible: 3
    }

    output {
        String regions_string = read_string("regions.txt")
    }
}

task SliceBam {
    input {
        String bam
        File bai
        String regions_args
        String sample_name
        Boolean is_terra = true
    }

    command <<<
        set -e

        apt-get update -qq && apt-get install -y -qq curl ca-certificates > /dev/null 2>&1

        if [ "~{is_terra}" = "true" ]; then
            export GCS_OAUTH_TOKEN=$(curl -s -H "Metadata-Flavor: Google" \
                http://metadata.google.internal/computeMetadata/v1/instance/service-accounts/default/token \
                | grep -o '"access_token":"[^"]*"' | cut -d'"' -f4)
            echo "GCS token acquired (${#GCS_OAUTH_TOKEN} chars)"
        fi

        echo "Slicing ~{sample_name}: ~{bam}"
        samtools view -@ 2 -X -b "~{bam}" "~{bai}" ~{regions_args} \
            | samtools sort -m 1G -@ 2 -o "~{sample_name}_mini.bam"
        samtools index "~{sample_name}_mini.bam"
        echo "Done: $(ls -lh ~{sample_name}_mini.bam)"
    >>>

    runtime {
        docker: "staphb/samtools:1.19"
        memory: "8 GB"
        cpu: 2
        disks: "local-disk 30 HDD"
        preemptible: 3
    }

    output {
        File minibam = "~{sample_name}_mini.bam"
        File minibai = "~{sample_name}_mini.bam.bai"
    }
}

task CreateScreenshots {
    input {
        File variants_tsv
        File p_bam
        File p_bai
        File m_bam
        File m_bai
        File f_bam
        File f_bai
        File ref_fasta
        File ref_index
        String out_name
        Int pad
        String docker_image
        Boolean is_terra = true
    }

    command <<<
        set -e

        mkdir screenshots

        # Clean input
        tr -d '\r' < "~{variants_tsv}" > variants.clean.tsv

        # Build IGV batch script
        cat > igv_batch.txt <<EOF
new
genome ~{ref_fasta}
load ~{p_bam}
load ~{m_bam}
load ~{f_bam}
snapshotDirectory screenshots
viewaspairs
EOF

        while IFS=$'\t' read -r chrom pos rest || [ -n "$chrom" ]; do
            if [[ ! "$pos" =~ ^[0-9]+$ ]]; then continue; fi
            if [[ "$chrom" != "chr"* ]]; then chrom="chr$chrom"; fi
            start=$((pos - ~{pad}))
            end=$((pos + ~{pad}))
            echo "goto ${chrom}:${start}-${end}" >> igv_batch.txt
            echo "sort base" >> igv_batch.txt
            echo "collapse" >> igv_batch.txt
            echo "snapshot ${chrom}_${pos}_trio.png" >> igv_batch.txt
        done < "variants.clean.tsv"

        echo "exit" >> igv_batch.txt

        echo "Running IGV on $(wc -l < variants.clean.tsv) variants..."
        igv_headless.sh -b igv_batch.txt

        echo "Archiving $(ls screenshots/*.png 2>/dev/null | wc -l) screenshots..."
        tar -czf ~{out_name}.tar.gz screenshots/
    >>>

    runtime {
        docker: docker_image
        memory: "6 GB"
        cpu: 2
        disks: "local-disk 20 HDD"
        preemptible: 3
    }

    output {
        File tar_gz = "~{out_name}.tar.gz"
    }
}

version 1.0

## MosaicPostProcessPart1: Trio VCF filtering + pileup validation
##
## Equivalent of getAllV6() + addPileup() from fanconi-genome.R:
##   Task 1 (FilterVariants): normalize trio VCFs, trio de novo filter (bcftools isec -C),
##           canonical chr + DP>=10 + segdup exclusion, extract positions + VCF info
##   Task 2 (TrioPileupNIO): samtools mpileup on trio BAMs via NIO streaming, merge with VCF info → CSV
##
## Part 2 (local): pileup-based de novo validation, hard filters, repeat/centromere, gnomAD, VAF binning

workflow MosaicPostProcessPart1 {
    input {
        File proband_vcf
        File proband_vcf_index
        File mother_vcf
        File mother_vcf_index
        File father_vcf
        File father_vcf_index

        String proband_bam
        String proband_bai
        String mother_bam
        String mother_bai
        String father_bam
        String father_bai

        File ref_fasta
        File ref_index

        File segdup_bed
        File lcr_bed

        String sample_name = "proband"
    }

    call FilterVariants {
        input:
            proband_vcf = proband_vcf,
            proband_vcf_index = proband_vcf_index,
            mother_vcf = mother_vcf,
            mother_vcf_index = mother_vcf_index,
            father_vcf = father_vcf,
            father_vcf_index = father_vcf_index,
            ref_fasta = ref_fasta,
            ref_index = ref_index,
            segdup_bed = segdup_bed,
            lcr_bed = lcr_bed,
            sample_name = sample_name
    }

    call TrioPileupNIO {
        input:
            positions_bed = FilterVariants.positions_bed,
            vcf_info_tsv = FilterVariants.vcf_info_tsv,
            proband_bam = proband_bam,
            proband_bai = proband_bai,
            mother_bam = mother_bam,
            mother_bai = mother_bai,
            father_bam = father_bam,
            father_bai = father_bai,
            ref_fasta = ref_fasta,
            ref_index = ref_index,
            sample_name = sample_name
    }

    output {
        File part1_csv = TrioPileupNIO.part1_csv
        File filter_report = FilterVariants.filter_report
    }
}

task FilterVariants {
    input {
        File proband_vcf
        File proband_vcf_index
        File mother_vcf
        File mother_vcf_index
        File father_vcf
        File father_vcf_index
        File ref_fasta
        File ref_index
        File segdup_bed
        File lcr_bed

        String sample_name

        String memory = "4 GB"
        Int cpu = 1
        Int disk_size = 30
        Int preemptible_attempts = 3
    }

    command <<<
        set -euo pipefail

        echo "=== FilterVariants: ~{sample_name} ===" > filter_report.txt
        echo "Input proband VCF: ~{proband_vcf}" >> filter_report.txt
        echo "Input mother VCF: ~{mother_vcf}" >> filter_report.txt
        echo "Input father VCF: ~{father_vcf}" >> filter_report.txt

        # Count raw variants
        RAW=$(bcftools view -H ~{proband_vcf} | wc -l)
        echo "Raw proband variants: ${RAW}" >> filter_report.txt

        # 1. Normalize all three VCFs: decompose multiallelics + left-align indels
        #    Then remove spanning deletions (ALT="*") — these are redundant markers
        #    for deletions already represented in other records
        echo "Normalizing proband VCF..."
        bcftools norm -m -any -f ~{ref_fasta} ~{proband_vcf} | \
            bcftools view -e 'ALT="*"' -O z -o proband_norm.vcf.gz
        bcftools index -t proband_norm.vcf.gz

        echo "Normalizing mother VCF..."
        bcftools norm -m -any -f ~{ref_fasta} ~{mother_vcf} | \
            bcftools view -e 'ALT="*"' -O z -o mother_norm.vcf.gz
        bcftools index -t mother_norm.vcf.gz

        echo "Normalizing father VCF..."
        bcftools norm -m -any -f ~{ref_fasta} ~{father_vcf} | \
            bcftools view -e 'ALT="*"' -O z -o father_norm.vcf.gz
        bcftools index -t father_norm.vcf.gz

        NORM=$(bcftools view -H proband_norm.vcf.gz | wc -l)
        echo "After proband normalize (excl spanning dels): ${NORM}" >> filter_report.txt

        # 2. Trio de novo filter: keep proband-private variants (not in mother or father)
        #    bcftools isec -C: complement — records from first file not in any other file
        #    Compares by CHROM, POS, REF, ALT
        bcftools isec -C -w 1 \
            proband_norm.vcf.gz mother_norm.vcf.gz father_norm.vcf.gz \
            -O z -o denovo.vcf.gz
        bcftools index -t denovo.vcf.gz

        DENOVO=$(bcftools view -H denovo.vcf.gz | wc -l)
        echo "After trio de novo filter (isec -C): ${DENOVO}" >> filter_report.txt

        # 3. Canonical chromosomes + DP >= 10
        CANONICAL="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY"
        bcftools view \
            -r ${CANONICAL} \
            -i 'FORMAT/DP>=10' \
            -O z \
            -o filtered.vcf.gz \
            denovo.vcf.gz
        bcftools index -t filtered.vcf.gz

        FILT=$(bcftools view -H filtered.vcf.gz | wc -l)
        echo "After canonical chr + DP>=10: ${FILT}" >> filter_report.txt

        # 4. Exclude segmental duplications
        bcftools view \
            -T ^~{segdup_bed} \
            -O z \
            -o nosegdup.vcf.gz \
            filtered.vcf.gz
        bcftools index -t nosegdup.vcf.gz

        NOSEGDUP=$(bcftools view -H nosegdup.vcf.gz | wc -l)
        echo "After segdup exclusion: ${NOSEGDUP}" >> filter_report.txt

        # 5. Exclude Low Complexity Regions (LCR)
        #    bcbio excluded LCR from HC calling intervals — our HC did not,
        #    producing ~730K false short indels in homopolymer/STR regions.
        #    Using Heng Li's LCR-hs38.bed.gz (1.75M regions, 228 Mb, 7.4% of genome)
        bcftools view \
            -T ^~{lcr_bed} \
            -O z \
            -o nolcr.vcf.gz \
            nosegdup.vcf.gz
        bcftools index -t nolcr.vcf.gz

        NOLCR=$(bcftools view -H nolcr.vcf.gz | wc -l)
        echo "After LCR exclusion: ${NOLCR}" >> filter_report.txt

        # 6. Extract positions BED for pileup
        bcftools query -f '%CHROM\t%POS0\t%END\n' nolcr.vcf.gz \
            | sort -k1,1 -k2,2n | uniq > positions.bed

        N_POS=$(wc -l < positions.bed)
        echo "Unique positions for pileup: ${N_POS}" >> filter_report.txt

        # 7. Extract VCF info TSV
        #    Missing values output as '.' by bcftools
        bcftools query \
            -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%INFO/MQ\t%INFO/QD\t%INFO/FS\t%INFO/ReadPosRankSum\t%INFO/BaseQRankSum\t%INFO/MQRankSum\t[%GT]\t[%AD]\t[%DP]\n' \
            nolcr.vcf.gz \
            > vcf_info.tsv

        echo "VCF info rows: $(wc -l < vcf_info.tsv)" >> filter_report.txt
        cat filter_report.txt
    >>>

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
        File positions_bed = "positions.bed"
        File vcf_info_tsv = "vcf_info.tsv"
        File filter_report = "filter_report.txt"
    }
}

task TrioPileupNIO {
    input {
        File positions_bed
        File vcf_info_tsv

        String proband_bam
        String proband_bai
        String mother_bam
        String mother_bai
        String father_bam
        String father_bai

        File ref_fasta
        File ref_index

        String sample_name

        String memory = "16 GB"
        Int cpu = 16
        Int disk_size = 30
        Int preemptible_attempts = 3
    }

    command <<<
        set -euo pipefail

        N_POS=$(wc -l < ~{positions_bed})
        echo "Positions to pileup: ${N_POS}"

        # Convert BED to GATK interval list format (chr:start-end, 1-based)
        awk '{print $1":"$2+1"-"$3}' ~{positions_bed} > positions.intervals
        echo "Intervals for PrintReads: $(wc -l < positions.intervals)"

        # Step 1: Extract BAM slices using gatk PrintReads (supports gs:// natively on Terra)
        # samtools cannot access gs:// URLs directly — GATK uses HTSJDK with built-in GCS support
        echo "Extracting proband BAM slice..."
        gatk PrintReads \
            -I ~{proband_bam} \
            --read-index ~{proband_bai} \
            -L positions.intervals \
            -O proband_slice.bam \
            --disable-read-filter NotDuplicateReadFilter \
            2>&1 | tail -5
        samtools index proband_slice.bam

        echo "Extracting mother BAM slice..."
        gatk PrintReads \
            -I ~{mother_bam} \
            --read-index ~{mother_bai} \
            -L positions.intervals \
            -O mother_slice.bam \
            --disable-read-filter NotDuplicateReadFilter \
            2>&1 | tail -5
        samtools index mother_slice.bam

        echo "Extracting father BAM slice..."
        gatk PrintReads \
            -I ~{father_bam} \
            --read-index ~{father_bai} \
            -L positions.intervals \
            -O father_slice.bam \
            --disable-read-filter NotDuplicateReadFilter \
            2>&1 | tail -5
        samtools index father_slice.bam

        echo "BAM slices ready:"
        ls -lh proband_slice.bam mother_slice.bam father_slice.bam

        # Step 2: Parallel per-position pileup on LOCAL BAM slices + merge with VCF info → CSV
        echo "Running parallel per-position pileup (${N_POS} positions, ~{cpu} workers)..."

        python3 <<'PYEOF'
import subprocess
import csv
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed

REF = "~{ref_fasta}"
BAMS = ["proband_slice.bam", "mother_slice.bam", "father_slice.bam"]
NCPU = int("~{cpu}")

def count_bases(pileup_str, ref_base):
    """Count bases from samtools mpileup pileup string."""
    counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'DEL': 0, 'INS': 0}
    i = 0
    s = pileup_str.upper()
    ref_upper = ref_base.upper()
    while i < len(s):
        c = s[i]
        if c in '.,':
            counts[ref_upper] = counts.get(ref_upper, 0) + 1
            i += 1
        elif c in 'ACGT':
            counts[c] += 1
            i += 1
        elif c in '*#':
            counts['DEL'] += 1
            i += 1
        elif c == '^':
            i += 2  # skip start-of-read + mapq char
        elif c == '$':
            i += 1
        elif c in '+-':
            if c == '+':
                counts['INS'] += 1
            else:
                counts['DEL'] += 1
            i += 1
            num_str = ''
            while i < len(s) and s[i].isdigit():
                num_str += s[i]
                i += 1
            if num_str:
                i += int(num_str)
        else:
            i += 1
    return counts

def pileup_one_position(chrom, pos):
    """Run samtools mpileup on a single position for all 3 local BAMs."""
    region = f"{chrom}:{pos}-{pos}"
    cmd = ["samtools", "mpileup", "-r", region, "-q", "20", "-Q", "10",
           "-d", "500", "-f", REF] + BAMS
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
    except subprocess.TimeoutExpired:
        sys.stderr.write(f"mpileup TIMEOUT at {region} (120s) — skipping\n")
        global n_timeouts
        n_timeouts += 1
        return chrom, pos, [0] * 21  # 3 samples × 7 fields
    if result.returncode != 0:
        sys.stderr.write(f"mpileup error at {region}: {result.stderr[:200]}\n")
        return chrom, pos, [0] * 21  # 3 samples × 7 fields

    sample_counts = []
    line = result.stdout.strip()
    if not line:
        return chrom, pos, [0] * 21

    parts = line.split('\t')
    ref = parts[2] if len(parts) > 2 else 'N'
    for sample_idx in range(3):
        offset = 3 + sample_idx * 3
        if offset >= len(parts):
            sample_counts.extend([0, 0, 0, 0, 0, 0, 0])
            continue
        depth = int(parts[offset])
        pileup_str = parts[offset + 1] if depth > 0 and offset + 1 < len(parts) else ''
        counts = count_bases(pileup_str, ref)
        sample_counts.extend([
            counts['A'], counts['C'], counts['G'], counts['T'],
            counts['DEL'], counts['INS'], depth
        ])
    return chrom, pos, sample_counts

# 1. Read positions from BED
positions = []
with open("~{positions_bed}") as f:
    for line in f:
        parts = line.strip().split('\t')
        chrom = parts[0]
        pos = str(int(parts[1]) + 1)  # BED 0-based → 1-based
        positions.append((chrom, pos))

print(f"Loaded {len(positions)} positions, using {NCPU} workers")

# 2. Parallel pileup on local BAM slices
pileup_data = {}
done = 0
n_errors = 0
n_timeouts = 0
with ThreadPoolExecutor(max_workers=NCPU) as executor:
    futures = {executor.submit(pileup_one_position, ch, p): (ch, p) for ch, p in positions}
    for future in as_completed(futures):
        chrom, pos, counts = future.result()
        pileup_data[f"{chrom}:{pos}"] = counts
        if all(c == 0 for c in counts):
            n_errors += 1
        done += 1
        if done % 1000 == 0:
            print(f"  Pileup progress: {done}/{len(positions)}")

print(f"Pileup done: {len(pileup_data)} positions ({n_errors} zero-count, {n_timeouts} timeouts)")

# 3. Read VCF info and merge with pileup → CSV
vcf_header = [
    'CHROM', 'POS', 'REF', 'ALT', 'QUAL',
    'MQ', 'QD', 'FS', 'ReadPosRankSum', 'BaseQRankSum', 'MQRankSum',
    'GT', 'AD', 'DP'
]

pileup_cols = []
for who in ['proband', 'mother', 'father']:
    for base in ['A', 'C', 'G', 'T', 'DEL', 'INS', 'depth']:
        pileup_cols.append(f"{who}_{base}")

out_header = vcf_header + pileup_cols
empty_pileup = [0] * len(pileup_cols)

n_matched = 0
n_missing = 0

with open("~{vcf_info_tsv}") as fin, \
     open('output.csv', 'w', newline='') as fout:
    writer = csv.writer(fout)
    writer.writerow(out_header)

    for line in fin:
        parts = line.strip().split('\t')
        if len(parts) < 14:
            continue
        chrom, pos = parts[0], parts[1]
        key = f"{chrom}:{pos}"
        pileup = pileup_data.get(key, empty_pileup)
        if key in pileup_data:
            n_matched += 1
        else:
            n_missing += 1
        writer.writerow(parts + [str(x) for x in pileup])

print(f"Matched: {n_matched}, Missing pileup: {n_missing}")
print(f"Total output rows: {n_matched + n_missing}")

PYEOF

        mv output.csv ~{sample_name}.part1.csv
        echo "Output: ~{sample_name}.part1.csv ($(wc -l < ~{sample_name}.part1.csv) lines incl header)"
    >>>

    runtime {
        docker: "broadinstitute/gatk:latest"
        memory: memory
        cpu: cpu
        disks: "local-disk " + disk_size + " SSD"
        preemptible: preemptible_attempts
        maxRetries: 3
        bootDiskSizeGb: 20
    }

    output {
        File part1_csv = "~{sample_name}.part1.csv"
    }
}

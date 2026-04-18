version 1.0

workflow TrioPileup {
    input {
        File proband_vcf
        File proband_vcf_index

        File proband_bam
        File proband_bai
        File mother_bam
        File mother_bai
        File father_bam
        File father_bai

        File ref_fasta
        File ref_index

        String sample_name = "proband"
        String memory = "8 GB"
        Int cpu = 2
        Int disk_size = 100
    }

    call NormalizeAndPileup {
        input:
            proband_vcf = proband_vcf,
            proband_vcf_index = proband_vcf_index,
            proband_bam = proband_bam,
            proband_bai = proband_bai,
            mother_bam = mother_bam,
            mother_bai = mother_bai,
            father_bam = father_bam,
            father_bai = father_bai,
            ref_fasta = ref_fasta,
            ref_index = ref_index,
            sample_name = sample_name,
            memory = memory,
            cpu = cpu,
            disk_size = disk_size
    }

    output {
        File normalized_vcf = NormalizeAndPileup.normalized_vcf
        File normalized_vcf_index = NormalizeAndPileup.normalized_vcf_index
        File pileup_tsv = NormalizeAndPileup.pileup_tsv
    }
}

task NormalizeAndPileup {
    input {
        File proband_vcf
        File proband_vcf_index

        File proband_bam
        File proband_bai
        File mother_bam
        File mother_bai
        File father_bam
        File father_bai

        File ref_fasta
        File ref_index

        String sample_name
        String memory
        Int cpu
        Int disk_size
    }

    command <<<
        set -euo pipefail

        # 1. Normalize VCF: decompose multiallelics + left-align indels
        bcftools norm \
            -m -any \
            -f ~{ref_fasta} \
            -O z \
            -o ~{sample_name}.norm.vcf.gz \
            ~{proband_vcf}
        bcftools index -t ~{sample_name}.norm.vcf.gz

        # 2. Extract variant positions for pileup (chrom:pos)
        bcftools query -f '%CHROM\t%POS\n' ~{sample_name}.norm.vcf.gz \
            | sort -k1,1 -k2,2n | uniq > positions.tsv

        N_VARIANTS=$(wc -l < positions.tsv)
        echo "Extracted ${N_VARIANTS} variant positions for pileup"

        # 3. Run samtools mpileup on trio at variant positions
        #    -l positions.tsv : restrict to positions
        #    -q 20 : min mapping quality
        #    -Q 10 : min base quality
        #    -d 1000 : max depth (for high-coverage data)
        #    -a : output all positions (even zero coverage)
        samtools mpileup \
            -l positions.tsv \
            -q 20 -Q 10 -d 1000 -a \
            -f ~{ref_fasta} \
            ~{proband_bam} ~{mother_bam} ~{father_bam} \
            > trio_pileup.raw.txt

        echo "Raw pileup lines: $(wc -l < trio_pileup.raw.txt)"

        # 4. Parse pileup into structured TSV
        #    mpileup output: CHROM POS REF [depth bases quals] x 3 samples
        python3 <<'PYEOF'
import sys
import re

def count_bases(pileup_str, ref_base):
    """Count bases from samtools mpileup pileup string."""
    counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'DEL': 0, 'INS': 0}
    i = 0
    s = pileup_str.upper()
    ref_upper = ref_base.upper()
    while i < len(s):
        c = s[i]
        if c == '.' or c == ',':
            # Reference base (. forward, , reverse - but we uppercased)
            counts[ref_upper] = counts.get(ref_upper, 0) + 1
            i += 1
        elif c in 'ACGT':
            counts[c] += 1
            i += 1
        elif c == '*' or c == '#':
            counts['DEL'] += 1
            i += 1
        elif c == '^':
            # Start of read + mapping quality char - skip both
            i += 2
        elif c == '$':
            # End of read
            i += 1
        elif c == '+' or c == '-':
            if c == '+':
                counts['INS'] += 1
            elif c == '-':
                counts['DEL'] += 1
            # Read the number following +/- and skip that many bases
            i += 1
            num_str = ''
            while i < len(s) and s[i].isdigit():
                num_str += s[i]
                i += 1
            if num_str:
                skip = int(num_str)
                i += skip
        else:
            i += 1
    total = sum(counts.values())
    return counts, total

header = ['CHROM', 'POS', 'REF',
          'proband_A', 'proband_C', 'proband_G', 'proband_T', 'proband_DEL', 'proband_INS', 'proband_depth',
          'mother_A', 'mother_C', 'mother_G', 'mother_T', 'mother_DEL', 'mother_INS', 'mother_depth',
          'father_A', 'father_C', 'father_G', 'father_T', 'father_DEL', 'father_INS', 'father_depth']

with open('trio_pileup.raw.txt') as fin, open('pileup_parsed.tsv', 'w') as fout:
    fout.write('\t'.join(header) + '\n')
    for line in fin:
        parts = line.strip().split('\t')
        if len(parts) < 12:
            continue
        chrom, pos, ref = parts[0], parts[1], parts[2]

        results = [chrom, pos, ref]
        for sample_idx in range(3):
            offset = 3 + sample_idx * 3
            depth = int(parts[offset])
            pileup_str = parts[offset + 1] if depth > 0 else ''
            counts, total = count_bases(pileup_str, ref)
            results.extend([
                str(counts['A']), str(counts['C']), str(counts['G']), str(counts['T']),
                str(counts['DEL']), str(counts['INS']), str(depth)
            ])
        fout.write('\t'.join(results) + '\n')

PYEOF

        # 5. Join pileup with VCF variant info (CHROM, POS, REF, ALT, QUAL, INFO fields)
        bcftools query \
            -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%INFO/MQ\t%INFO/QD\t%INFO/FS\t%INFO/ReadPosRankSum\t%INFO/BaseQRankSum\t%INFO/MQRankSum\t[%GT]\t[%AD]\t[%DP]\n' \
            ~{sample_name}.norm.vcf.gz \
            > vcf_info.tsv

        python3 <<'PYEOF2'
# Join VCF info with pileup data on CHROM:POS
vcf_info = {}
with open('vcf_info.tsv') as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) < 14:
            continue
        chrom, pos = parts[0], parts[1]
        key = f"{chrom}:{pos}"
        # May have multiple ALTs at same position after norm - keep all
        if key not in vcf_info:
            vcf_info[key] = []
        vcf_info[key].append(parts)

pileup_data = {}
with open('pileup_parsed.tsv') as f:
    header_line = f.readline().strip()
    pileup_cols = header_line.split('\t')[3:]  # skip CHROM, POS, REF
    for line in f:
        parts = line.strip().split('\t')
        chrom, pos = parts[0], parts[1]
        key = f"{chrom}:{pos}"
        pileup_data[key] = parts[3:]

out_header = [
    'CHROM', 'POS', 'REF', 'ALT', 'QUAL',
    'MQ', 'QD', 'FS', 'ReadPosRankSum', 'BaseQRankSum', 'MQRankSum',
    'GT', 'AD', 'DP'
] + pileup_cols

with open('output.tsv', 'w') as fout:
    fout.write('\t'.join(out_header) + '\n')
    for key, vcf_rows in vcf_info.items():
        pileup = pileup_data.get(key, ['0'] * len(pileup_cols))
        for row in vcf_rows:
            out_row = row + pileup
            fout.write('\t'.join(out_row) + '\n')

PYEOF2

        mv output.tsv ~{sample_name}.pileup.tsv
        echo "Final output rows: $(wc -l < ~{sample_name}.pileup.tsv)"
    >>>

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
        File normalized_vcf = "~{sample_name}.norm.vcf.gz"
        File normalized_vcf_index = "~{sample_name}.norm.vcf.gz.tbi"
        File pileup_tsv = "~{sample_name}.pileup.tsv"
    }
}

#!/usr/bin/env python3
"""
Cross-caller comparison: HC Part2 vs Mutect2 Part2 vs Old Pipeline (bcbio HC).

Creates:
  1. Master CSV with all variants and flags per caller/filter
  2. Intersection tables (mosaic SNV, het SNV)
  3. Venn diagrams (matplotlib_venn)

Usage:
  python scripts/cross_caller_comparison.py
"""

import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from _common import (
    OUTPUTS_DIR, PROJECT_DIR, TRIOS,
    BH_SIG_STRONG, BH_SIG_MODERATE, classify_vaf,
)

try:
    from matplotlib_venn import venn3
    HAS_VENN = True
except ImportError:
    HAS_VENN = False
    print("WARNING: matplotlib_venn not installed. pip install matplotlib-venn")

# --- Paths ---
PART2_DIR = OUTPUTS_DIR / "part2"
OLD_PIPELINE_DIR = PROJECT_DIR / "data" / "metadata"
OUTPUT_DIR = OUTPUTS_DIR / "cross_caller"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def _safe_float(val):
    try:
        return float(val)
    except (ValueError, TypeError):
        return None


def load_hc_part2(trio):
    """Load HC Part2 final variants with all columns."""
    f = PART2_DIR / f"{trio}_part2_final.tsv"
    variants = {}
    with open(f) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            key = f"{row['CHROM']}:{row['POS']}:{row['REF']}:{row['ALT']}"
            vaf = float(row['VAF'])
            vtype = "SNV" if (len(row['REF']) == 1 and len(row['ALT']) == 1) else "INDEL"
            d = {
                'CHROM': row['CHROM'], 'POS': row['POS'],
                'REF': row['REF'], 'ALT': row['ALT'],
                'VAF_hc': vaf, 'vtype': vtype, 'vaf_class_hc': classify_vaf(vaf),
                'deepSNV_BH_hc': _safe_float(row.get('deepSNV_pval_BH', '')),
                '_raw_hc': row,  # store full row
            }
            variants[key] = d
    return variants


def load_mutect2_part2(trio):
    """Load Mutect2 Part2 final variants with all columns."""
    f = PART2_DIR / f"{trio}_mutect2_part2_final.tsv"
    if not f.exists():
        return {}
    variants = {}
    with open(f) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            key = f"{row['CHROM']}:{row['POS']}:{row['REF']}:{row['ALT']}"
            vaf = float(row['VAF'])
            vtype = "SNV" if (len(row['REF']) == 1 and len(row['ALT']) == 1) else "INDEL"
            d = {
                'CHROM': row['CHROM'], 'POS': row['POS'],
                'REF': row['REF'], 'ALT': row['ALT'],
                'VAF_m2': vaf, 'vtype': vtype, 'vaf_class_m2': classify_vaf(vaf),
                'deepSNV_BH_m2': _safe_float(row.get('deepSNV_pval_BH', '')),
                '_raw_m2': row,  # store full row
            }
            variants[key] = d
    return variants


def load_old_pipeline(trio):
    """Load old pipeline (bcbio HC) variants."""
    f = OLD_PIPELINE_DIR / f"all-vars-{trio}.csv"
    variants = {}
    with open(f) as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            chrom = row['Chr']
            if not chrom.startswith('chr'):
                chrom = f"chr{chrom}"
            pos = row['POS']
            ref = row['Ref']
            alt = row['Alt']
            key = f"{chrom}:{pos}:{ref}:{alt}"
            vaf = float(row['VAF'])
            vtype = "SNV" if (len(ref) == 1 and len(alt) == 1) else "INDEL"
            variants[key] = {
                'CHROM': chrom, 'POS': pos, 'REF': ref, 'ALT': alt,
                'VAF_old': vaf, 'vtype': vtype, 'vaf_class_old': classify_vaf(vaf),
            }
    return variants


def make_master_csv(trio, hc, m2, old):
    """Create master CSV with all variants, flags, and detailed annotations."""
    all_keys = set(hc.keys()) | set(m2.keys()) | set(old.keys())

    # HC detail columns to include
    HC_DETAIL_COLS = [
        'QUAL', 'MQ', 'QD', 'FS', 'ReadPosRankSum', 'BaseQRankSum', 'MQRankSum',
        'GT', 'AD', 'DP',
        'proband_A', 'proband_C', 'proband_G', 'proband_T', 'proband_DEL', 'proband_INS', 'proband_depth',
        'mother_A', 'mother_C', 'mother_G', 'mother_T', 'mother_DEL', 'mother_INS', 'mother_depth',
        'father_A', 'father_C', 'father_G', 'father_T', 'father_DEL', 'father_INS', 'father_depth',
        'AD_alt', 'AD_ref', 'gnomADg_AF',
    ]
    # Mutect2-specific detail columns
    M2_DETAIL_COLS = [
        'QUAL', 'TLOD', 'GERMQ', 'SEQQ', 'STRANDQ', 'MBQ', 'MPOS',
        'GT', 'AD', 'DP', 'AF',
        'proband_A', 'proband_C', 'proband_G', 'proband_T', 'proband_DEL', 'proband_INS', 'proband_depth',
        'mother_A', 'mother_C', 'mother_G', 'mother_T', 'mother_DEL', 'mother_INS', 'mother_depth',
        'father_A', 'father_C', 'father_G', 'father_T', 'father_DEL', 'father_INS', 'father_depth',
        'AD_alt', 'AD_ref', 'gnomADg_AF',
    ]

    rows = []
    for key in sorted(all_keys):
        parts = key.split(':')
        chrom, pos, ref, alt = parts[0], parts[1], parts[2], parts[3]
        vtype = "SNV" if (len(ref) == 1 and len(alt) == 1) else "INDEL"

        in_hc = key in hc
        in_m2 = key in m2
        in_old = key in old

        # Only include variants that are mosaic or het in at least one caller
        classes = []
        if in_hc: classes.append(hc[key].get('vaf_class_hc', ''))
        if in_m2: classes.append(m2[key].get('vaf_class_m2', ''))
        if in_old: classes.append(old[key].get('vaf_class_old', ''))
        if not any(c in ('mosaic', 'het') for c in classes):
            continue

        vaf_hc = hc[key]['VAF_hc'] if in_hc else None
        vaf_m2 = m2[key]['VAF_m2'] if in_m2 else None
        vaf_old = old[key]['VAF_old'] if in_old else None

        vaf_class_hc = hc[key].get('vaf_class_hc', '') if in_hc else ''
        vaf_class_m2 = m2[key].get('vaf_class_m2', '') if in_m2 else ''
        vaf_class_old = old[key].get('vaf_class_old', '') if in_old else ''
        consensus_class = vaf_class_hc or vaf_class_m2 or vaf_class_old

        bh_hc = hc[key].get('deepSNV_BH_hc') if in_hc else None
        bh_m2 = m2[key].get('deepSNV_BH_m2') if in_m2 else None

        # Union column: HC ∪ Mutect2 (in at least one of the new callers)
        in_union = in_hc or in_m2

        row = {
            'CHROM': chrom, 'POS': pos, 'REF': ref, 'ALT': alt,
            'variant_type': vtype,
            'in_HC': in_hc, 'in_Mutect2': in_m2, 'in_HC_or_M2': in_union, 'in_OldPipeline': in_old,
            'VAF_HC': f"{vaf_hc:.4f}" if vaf_hc else '',
            'VAF_Mutect2': f"{vaf_m2:.4f}" if vaf_m2 else '',
            'VAF_OldPipeline': f"{vaf_old:.4f}" if vaf_old else '',
            'class_HC': vaf_class_hc,
            'class_Mutect2': vaf_class_m2,
            'class_OldPipeline': vaf_class_old,
            'consensus_class': consensus_class,
            'deepSNV_BH_HC': f"{bh_hc:.2e}" if bh_hc is not None else '',
            'deepSNV_BH_Mutect2': f"{bh_m2:.2e}" if bh_m2 is not None else '',
            'deepSNV_BH_HC_sig001': bh_hc is not None and bh_hc < BH_SIG_STRONG if in_hc else '',
            'deepSNV_BH_HC_sig005': bh_hc is not None and bh_hc < BH_SIG_MODERATE if in_hc else '',
            'deepSNV_BH_M2_sig001': bh_m2 is not None and bh_m2 < BH_SIG_STRONG if in_m2 else '',
            'deepSNV_BH_M2_sig005': bh_m2 is not None and bh_m2 < BH_SIG_MODERATE if in_m2 else '',
        }

        # Add HC detail columns (prefixed with HC_)
        raw_hc = hc[key].get('_raw_hc', {}) if in_hc else {}
        for col in HC_DETAIL_COLS:
            row[f'HC_{col}'] = raw_hc.get(col, '')

        # Add Mutect2 detail columns (prefixed with M2_)
        raw_m2 = m2[key].get('_raw_m2', {}) if in_m2 else {}
        for col in M2_DETAIL_COLS:
            row[f'M2_{col}'] = raw_m2.get(col, '')

        rows.append(row)

    outpath = OUTPUT_DIR / f"{trio}_master.csv"
    fieldnames = list(rows[0].keys()) if rows else []
    with open(outpath, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    return rows, outpath


def compute_intersections(hc, m2, old, vtype_filter="SNV", class_filter="mosaic"):
    """Compute set intersections for a given variant type and VAF class."""
    def get_keys(d, caller_class_field):
        return {k for k, v in d.items()
                if v['vtype'] == vtype_filter and v.get(caller_class_field, '') == class_filter}

    hc_set = get_keys(hc, 'vaf_class_hc')
    m2_set = get_keys(m2, 'vaf_class_m2')
    old_set = get_keys(old, 'vaf_class_old')

    return hc_set, m2_set, old_set


def venn_diagram(hc_set, m2_set, old_set, title, outpath):
    """Create Venn diagram."""
    if not HAS_VENN:
        return

    fig, ax = plt.subplots(figsize=(8, 6))
    v = venn3([hc_set, m2_set, old_set],
              set_labels=('HC (new)', 'Mutect2', 'Old HC (bcbio)'),
              ax=ax)
    ax.set_title(title, fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.close()


def deepsnv_intersection(hc, m2, vtype_filter="SNV", class_filter="mosaic"):
    """Get deepSNV-significant variants in shared HC∩M2 set."""
    hc_keys = {k for k, v in hc.items() if v['vtype'] == vtype_filter and v.get('vaf_class_hc', '') == class_filter}
    m2_keys = {k for k, v in m2.items() if v['vtype'] == vtype_filter and v.get('vaf_class_m2', '') == class_filter}
    shared = hc_keys & m2_keys

    sig_hc_005 = sum(1 for k in shared if hc[k].get('deepSNV_BH_hc') is not None and hc[k]['deepSNV_BH_hc'] < BH_SIG_MODERATE)
    sig_hc_001 = sum(1 for k in shared if hc[k].get('deepSNV_BH_hc') is not None and hc[k]['deepSNV_BH_hc'] < BH_SIG_STRONG)
    sig_m2_005 = sum(1 for k in shared if m2[k].get('deepSNV_BH_m2') is not None and m2[k]['deepSNV_BH_m2'] < BH_SIG_MODERATE)
    sig_m2_001 = sum(1 for k in shared if m2[k].get('deepSNV_BH_m2') is not None and m2[k]['deepSNV_BH_m2'] < BH_SIG_STRONG)

    return len(shared), sig_hc_001, sig_hc_005, sig_m2_001, sig_m2_005


def main():
    summary_lines = []

    for trio in TRIOS:
        print(f"\n{'='*60}")
        print(f"  {trio.upper()}")
        print(f"{'='*60}")

        hc = load_hc_part2(trio)
        m2 = load_mutect2_part2(trio)
        old = load_old_pipeline(trio)

        # Master CSV
        rows, master_path = make_master_csv(trio, hc, m2, old)
        print(f"Master CSV: {master_path} ({len(rows)} variants)")

        for class_filter in ["mosaic", "het"]:
            for vtype in ["SNV", "INDEL"]:
                hc_set, m2_set, old_set = compute_intersections(hc, m2, old, vtype, class_filter)

                if not (hc_set or m2_set or old_set):
                    continue

                shared_all = hc_set & m2_set & old_set
                hc_m2 = hc_set & m2_set
                hc_old = hc_set & old_set
                m2_old = m2_set & old_set
                hc_only = hc_set - m2_set - old_set
                m2_only = m2_set - hc_set - old_set
                old_only = old_set - hc_set - m2_set

                label = f"{class_filter} {vtype}"
                print(f"\n  --- {label} ---")
                print(f"  HC={len(hc_set)}, M2={len(m2_set)}, Old={len(old_set)}")
                print(f"  HC∩M2={len(hc_m2)}, HC∩Old={len(hc_old)}, M2∩Old={len(m2_old)}, All3={len(shared_all)}")
                print(f"  HC-only={len(hc_only)}, M2-only={len(m2_only)}, Old-only={len(old_only)}")

                # deepSNV for shared HC∩M2
                if class_filter in ["mosaic", "het"] and vtype == "SNV":
                    n_shared, sig001_hc, sig005_hc, sig001_m2, sig005_m2 = deepsnv_intersection(
                        hc, m2, vtype, class_filter)
                    if n_shared > 0:
                        print(f"  Shared HC∩M2 deepSNV: HC BH<0.05={sig005_hc}/{n_shared}, HC BH<0.01={sig001_hc}/{n_shared}")
                        print(f"                        M2 BH<0.05={sig005_m2}/{n_shared}, M2 BH<0.01={sig001_m2}/{n_shared}")

                # Venn diagram
                if hc_set or m2_set or old_set:
                    venn_path = OUTPUT_DIR / f"{trio}_{class_filter}_{vtype}_venn.png"
                    venn_diagram(hc_set, m2_set, old_set,
                                f"{trio} — {label}", venn_path)

        # Summary table row
        hc_mosaic_snv = {k for k, v in hc.items() if v['vtype'] == 'SNV' and v.get('vaf_class_hc') == 'mosaic'}
        m2_mosaic_snv = {k for k, v in m2.items() if v['vtype'] == 'SNV' and v.get('vaf_class_m2') == 'mosaic'}
        old_mosaic_snv = {k for k, v in old.items() if v['vtype'] == 'SNV' and v.get('vaf_class_old') == 'mosaic'}
        hc_het_snv = {k for k, v in hc.items() if v['vtype'] == 'SNV' and v.get('vaf_class_hc') == 'het'}
        m2_het_snv = {k for k, v in m2.items() if v['vtype'] == 'SNV' and v.get('vaf_class_m2') == 'het'}
        old_het_snv = {k for k, v in old.items() if v['vtype'] == 'SNV' and v.get('vaf_class_old') == 'het'}

        shared_mosaic = hc_mosaic_snv & m2_mosaic_snv
        shared_het = hc_het_snv & m2_het_snv

        # deepSNV on shared mosaic
        n_sh, sig001_h, sig005_h, sig001_m, sig005_m = deepsnv_intersection(hc, m2, "SNV", "mosaic")

        summary_lines.append({
            'trio': trio,
            'hc_mosaic': len(hc_mosaic_snv),
            'm2_mosaic': len(m2_mosaic_snv),
            'old_mosaic': len(old_mosaic_snv),
            'shared_hc_m2_mosaic': len(shared_mosaic),
            'shared_all3_mosaic': len(hc_mosaic_snv & m2_mosaic_snv & old_mosaic_snv),
            'hc_het': len(hc_het_snv),
            'm2_het': len(m2_het_snv),
            'old_het': len(old_het_snv),
            'shared_hc_m2_het': len(shared_het),
            'shared_all3_het': len(hc_het_snv & m2_het_snv & old_het_snv),
            'shared_mosaic_deepSNV_hc_005': sig005_h,
            'shared_mosaic_deepSNV_m2_005': sig005_m,
        })

    # Print summary tables
    print(f"\n{'='*80}")
    print("CROSS-CALLER SUMMARY — MOSAIC SNV")
    print(f"{'='*80}")
    print(f"{'Sample':<10} {'HC':>5} {'M2':>5} {'Old':>5} {'HC∩M2':>6} {'All3':>5} {'HC∩M2 sig(HC BH<.05)':>22} {'HC∩M2 sig(M2 BH<.05)':>22}")
    for s in summary_lines:
        sh = s['shared_hc_m2_mosaic']
        print(f"{s['trio']:<10} {s['hc_mosaic']:>5} {s['m2_mosaic']:>5} {s['old_mosaic']:>5} {sh:>6} {s['shared_all3_mosaic']:>5} {s['shared_mosaic_deepSNV_hc_005']:>3}/{sh:<3} {'':<14} {s['shared_mosaic_deepSNV_m2_005']:>3}/{sh:<3}")

    print(f"\n{'='*80}")
    print("CROSS-CALLER SUMMARY — HET DE NOVO SNV")
    print(f"{'='*80}")
    print(f"{'Sample':<10} {'HC':>5} {'M2':>5} {'Old':>5} {'HC∩M2':>6} {'All3':>5}")
    for s in summary_lines:
        print(f"{s['trio']:<10} {s['hc_het']:>5} {s['m2_het']:>5} {s['old_het']:>5} {s['shared_hc_m2_het']:>6} {s['shared_all3_het']:>5}")

    print(f"\nOutput directory: {OUTPUT_DIR}")
    print(f"Master CSVs: {OUTPUT_DIR}/{{trio}}_master.csv")
    if HAS_VENN:
        print(f"Venn diagrams: {OUTPUT_DIR}/{{trio}}_{{class}}_{{type}}_venn.png")


if __name__ == "__main__":
    main()

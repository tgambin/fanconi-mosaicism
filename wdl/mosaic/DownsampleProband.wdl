version 1.0

workflow DownsampleProband {
    input {
        File proband_bam
        File proband_bai
        String sample_name
        Array[String] target_coverages = ["150", "100", "50"]

        File? ref_fasta
        File? ref_index
        Boolean is_wgs = true
    }

    call Downsample {
        input:
            bam = proband_bam,
            bai = proband_bai,
            sample_name = sample_name,
            target_coverages = target_coverages,
            is_wgs = is_wgs
    }
    
    output {
        Array[File] downsampled_bams = Downsample.output_bams
        Array[File] downsampled_bais = Downsample.output_bais
        Array[File] coverage_reports = Downsample.output_coverage_reports
        File initial_coverage_report = Downsample.initial_coverage_report
    }
}

task Downsample {
    input {
        File bam
        File bai
        String sample_name
        Array[String] target_coverages
        Boolean is_wgs
    }


    # Estimate disk size to accommodate input BAM (300G) + outputs (150G + 100G + 50G ~ 300G) + temp
    # 1000G to be safe and allow for all outputs to be kept for WDL glob
    Int disk_size = 1000

    command <<<
        set -e
        
        # 0. Setup Environment: Install Python3, pip
        echo "Installing dependencies..."
        if ! command -v python3 &> /dev/null; then
            apt-get update && apt-get install -y python3 python3-pip curl
        fi
        
        # Function to calculate coverage stats
        calc_stats() {
            local input_bam="$1"
            local output_txt="$2"
            
            # Determine samtools depth flags based on WGS vs Slice
            # If is_wgs is TRUE, use -a (all positions including zero coverage)
            # If is_wgs is FALSE (slice), do NOT use -a (only covered positions)
            local depth_flags="-@ 7"
            if [ "~{is_wgs}" == "true" ]; then
                depth_flags="-@ 7 -a"
                echo "Running as WGS: using samtools depth -a"
            else
                echo "Running as Slice: omitting -a from samtools depth"
            fi
            
            echo "Calculating coverage for ${input_bam}..."
            samtools depth ${depth_flags} "${input_bam}" | \
            awk '{sum+=$3; sumsq+=$3*$3; count++} END { 
                if (count > 0) { 
                    print "Average = ", sum/count; 
                    print "Stdev = ", sqrt(sumsq/count - (sum/count)^2) 
                } else { 
                    print "Average = 0"; 
                    print "Stdev = 0" 
                } 
            }' > "${output_txt}"
        }

        # Function for downsampling, stats calculation
        process_level() {
            local target_cov="$1"   # Target coverage (e.g., 150.0)
            local baseline_cov="$2" # Current average coverage
            local suffix="$3"       # Output suffix (e.g., "150x")
            
            # Check if downsampling is needed (baseline > target) using python for float comparison
            local needed=$(python3 -c "print(1 if ${baseline_cov} > ${target_cov} else 0)")
            
            if [ "$needed" -eq 1 ]; then
                # Calculate fraction and seed
                # Fraction = target / baseline (capped at 1.0)
                # Seed = 42 + fractional_part
                python3 -c "
fraction = ${target_cov} / ${baseline_cov}
if fraction > 1.0:
    fraction = 1.0
# Ensure seed part is 42, fraction part is the fraction
# Output format: FRACTION SEED.FRACTION
seed_arg = 42 + fraction
print(f'{fraction} {seed_arg:.6f}')
" > params.txt
                
                local fraction=$(awk '{print $1}' params.txt)
                local seed_arg=$(awk '{print $2}' params.txt)
                local out_bam="~{sample_name}.${suffix}.bam"
                local out_cov="~{sample_name}.${suffix}_coverage.txt"
                
                echo "Downsampling ~{bam} to ${target_cov}x from ${baseline_cov}x (Fraction: ${fraction})"
                
                # 1. Downsample
                samtools view -@ 7 -b -s "$seed_arg" -o "$out_bam" ~{bam}
                samtools index -@ 7 "$out_bam"
                
                # 2. Calculate Coverage
                calc_stats "$out_bam" "$out_cov"
                
                # 3. Output files for collection
                # WDL output block will pick these up via glob
                
            else
                echo "Skipping downsampling to ${target_cov}x (Baseline ${baseline_cov}x is sufficient or lower)"
            fi
        }

        # --- Main Execution ---

        # 1. Calculate Initial Coverage
        initial_stats="~{sample_name}.initial_coverage.txt"
        calc_stats "~{bam}" "$initial_stats"
        
        # Parse average coverage
        avg_cov=$(grep "Average =" "$initial_stats" | awk '{print $3}')
        echo "Initial Average Coverage: $avg_cov"
        
        # 2. Perform Downsampling Levels
        TARGETS="~{sep=' ' target_coverages}"
        for t in $TARGETS; do
            process_level "${t}.0" "$avg_cov" "${t}x"
        done
        
        echo "All requested downsampling completed."
    >>>

    runtime {
        docker: "staphb/samtools:1.19"
        cpu: 8
        memory: "32 GB"
        disks: "local-disk " + disk_size + " SSD"
        preemptible: 0
    }
    
    output {
        File initial_coverage_report = "~{sample_name}.initial_coverage.txt"
        Array[File] output_bams = glob("*.bam")
        Array[File] output_bais = glob("*.bam.bai")
        Array[File] output_coverage_reports = glob("*_coverage.txt")
    }
}

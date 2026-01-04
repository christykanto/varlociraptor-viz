import pysam
import pandas as pd
import altair as alt
alt.data_transformers.disable_max_rows()

# -------------------------
# Helper
# -------------------------
def phred_to_prob(q):
    return 10 ** (-q / 10)

# -------------------------
# 1. Event probabilities
# -------------------------
def plot_event_probabilities(record):
    data = []
    for key, value in record.info.items():
        if key.startswith("PROB_"):
            prob = phred_to_prob(float(value))
            data.append({
                "event": key.replace("PROB_", ""),
                "probability": prob,
                "phred": float(value)
            })
    
    df = pd.DataFrame(data)
    print("Event probabilities:")
    print(df)
    
    # Use linear scale instead of log
    chart = (
        alt.Chart(df)
        .mark_bar()
        .encode(
            x=alt.X("event:N", title="Event"),
            y=alt.Y(
                "probability:Q",
                title="Probability"
                # No log scale - use linear
            ),
            color=alt.Color("event:N"),
            tooltip=["event", alt.Tooltip("probability:Q", format=".4f"), "phred"]
        )
        .properties(title="Event Probabilities (PHRED converted)", width=400, height=300)
    )
    return chart

# -------------------------
# 2. Allele frequency distribution
# -------------------------
def plot_allele_frequency_distribution(record, sample):
    sample_data = record.samples[sample]
    afd = sample_data["AFD"]          # FORMAT field
    mle_af = float(sample_data["AF"]) # FORMAT field
    
    rows = []
    for entry in afd:
        af, q = entry.split("=")
        rows.append({
            "allele_frequency": float(af),
            "probability": phred_to_prob(float(q)),
            "is_mle": abs(float(af) - mle_af) < 0.001  # Use tolerance for comparison
        })
    
    df = pd.DataFrame(rows)
    print("\nAFD data points:", len(df))
    
    chart = (
        alt.Chart(df)
        .mark_circle(size=100)
        .encode(
            x=alt.X("allele_frequency:Q", title="Allele Frequency"),
            y=alt.Y(
                "probability:Q",
                title="Probability",
                scale=alt.Scale(type="log")
            ),
            color=alt.condition(
                "datum.is_mle",
                alt.value("red"),
                alt.value("steelblue")
            ),
            tooltip=["allele_frequency", "probability", "is_mle"]
        )
        .properties(title="Allele Frequency Distribution", width=500, height=350)
    )
    return chart

# -------------------------
# 3. Observations (OBS)
# -------------------------
import re

def parse_obs_detailed(obs_list):
    """
    Parse OBS entries in detail.
    Format: CBDTASOPXI
    C = count (digits)
    B = 2-letter posterior odds code
    D = edit distance
    T = alignment type (s/p)
    A = alt locus (#/*./)
    S = strand (+/-/*)
    O = orientation (>/</*!/!)
    P = read position (^/*)
    X = softclip ($/.)
    I = indel (*/.)
    """
    parsed = []
    
    for obs in obs_list:
        # Match: digits + 2 letters + 8 single chars
        match = re.match(r'^(\d+)([ARar][A-Za-z])(.)(.)(.)(.)(.)(.)(.)(.)$', obs)
        if match:
            count, odds_code, edit_dist, align_type, alt_locus, strand, orient, read_pos, softclip, indel = match.groups()
            
            # Determine allele
            allele = "REF" if odds_code[0].upper() == 'R' else "ALT"
            
            parsed.append({
                "allele": allele,
                "count": int(count),
                "odds_code": odds_code,
                "edit_distance": 0 if edit_dist == '.' else int(edit_dist),
                "alignment_type": align_type,
                "alt_locus": alt_locus,
                "strand": strand,
                "orientation": orient,
                "read_position": read_pos,
                "softclip": softclip,
                "indel": indel
            })
    
    return pd.DataFrame(parsed)

def plot_observations(record, sample):
    obs_list = record.samples[sample]["OBS"]
    obs_df = parse_obs_detailed(obs_list)
    
    print(f"\nParsed {len(obs_df)} OBS entries")
    
    # Prepare data - each observation gets ALL components
    data_for_plot = []
    
    for idx, row in obs_df.iterrows():
        obs_id = f"obs_{idx}"
        
        # Each component type as separate entry but same observation
        components = [
            ("Posterior Odds", row["odds_code"]),
            ("Edit Distance", f"dist={row['edit_distance']}" if row['edit_distance'] > 0 else "."),
            ("Strand", row["strand"]),
            ("Orientation", row["orientation"]),
            ("Read Position", row["read_position"]),
            ("Softclip", row["softclip"]),
            ("Indel", row["indel"])
        ]
        
        for comp_name, comp_value in components:
            data_for_plot.append({
                "allele": row["allele"],
                "observation": obs_id,
                "component": comp_name,
                "category": comp_value,
                "value": row["count"]
            })
    
    plot_df = pd.DataFrame(data_for_plot)
    
    # Create chart with components side-by-side within each allele
    chart = (
        alt.Chart(plot_df)
        .mark_bar()
        .encode(
            x=alt.X("component:N", title="Component Type", axis=alt.Axis(labelAngle=-45)),
            y=alt.Y("value:Q", title="Count", stack="zero"),
            color=alt.Color("category:N", title="Category"),
            row=alt.Row("allele:N", title="Allele Type"),
            tooltip=["allele", "observation", "component", "category", "value"]
        )
        .properties(width=600, height=250)
    )
    
    return chart

# -------------------------
# Main
# -------------------------
def main():
    vcf_path = "/home/christy/varlociraptor-viz-repo/examples/example.vcf"
    sample_name = "sample1"
    
    vcf = pysam.VariantFile(vcf_path)
    record = next(iter(vcf))
    
    print(f"Processing variant at {record.chrom}:{record.pos}")
    print(f"REF: {record.ref}, ALT: {record.alts}")
    
    plot_event_probabilities(record).save("event_probs.html")
    plot_allele_frequency_distribution(record, sample_name).save("afd.html")
    plot_observations(record, sample_name).save("obs.html")
    
    print("\nPlots saved:")
    print(" - event_probs.html")
    print(" - afd.html")
    print(" - obs.html")

if __name__ == "__main__":
    main()

import altair as alt
import pandas as pd
import pysam
import re

def phred_to_prob(phred_value):
    """Convert PHRED score to probability"""
    # Handle tuples (take first value) or single values
    if isinstance(phred_value, (tuple, list)):
        phred_value = phred_value[0]
    return 10 ** (-phred_value / 10)

def visualize_event_probabilities(record):
    """
    Visualize event probabilities from INFO column (PROB_* fields)
    """
    # Extract all PROB_ fields from INFO
    prob_data = []
    for key, value in record.info.items():
        if key.startswith('PROB_'):
            event_name = key.replace('PROB_', '')
            probability = phred_to_prob(value)
            prob_data.append({'Event': event_name, 'Probability': probability})
    
    # Create DataFrame
    df = pd.DataFrame(prob_data)
    
    # Create bar plot with log scale
    chart = alt.Chart(df).mark_bar().encode(
        x=alt.X('Event:N', title='Event Type'),
        y=alt.Y('Probability:Q', scale=alt.Scale(type='log'), title='Probability (log scale)'),
        tooltip=['Event', 'Probability']
    ).properties(
        title='Event Probabilities',
        width=400,
        height=300
    )
    
    return chart

def visualize_allele_frequency_distribution(record, sample_name):
    """
    Visualize allele frequency distribution (AFD field)
    ML estimate is layered on top of distribution points
    """
    # Get sample data
    sample = record.samples[sample_name]
    
    # Parse AFD field (allele frequency distribution)
    afd = sample['AFD']
    af_ml = sample['AF']  # Maximum likelihood allele frequency
    
    # Handle if AF is a tuple
    if isinstance(af_ml, (tuple, list)):
        af_ml = af_ml[0]
    
    # Parse AFD entries: format is "freq1=phred1,freq2=phred2,..."
    afd_data = []
    
    # Handle if AFD is a tuple or list
    if isinstance(afd, (tuple, list)):
        afd_entries = afd
    else:
        afd_entries = [afd]
    
    for entry in afd_entries:
        # Each entry might be a string with multiple comma-separated values
        if isinstance(entry, str):
            parts = entry.split(',')
            for part in parts:
                if '=' in part:
                    freq, phred = part.split('=')
                    freq = float(freq)
                    prob = phred_to_prob(float(phred))
                    is_ml = (abs(freq - af_ml) < 0.001)
                    afd_data.append({
                        'Allele Frequency': freq,
                        'Probability': prob,
                        'Type': 'ML Estimate' if is_ml else 'Distribution'
                    })
    
    df = pd.DataFrame(afd_data)
    
    # Split data into distribution and ML estimate
    df_distribution = df[df['Type'] == 'Distribution']
    df_ml = df[df['Type'] == 'ML Estimate']
    
    # Create base layer for distribution points (blue)
    base_layer = alt.Chart(df_distribution).mark_circle(size=60, opacity=0.7).encode(
        x=alt.X('Allele Frequency:Q', title='Allele Frequency'),
        y=alt.Y('Probability:Q', axis=None),
        color=alt.value('blue'),
        tooltip=['Allele Frequency', 'Probability', 'Type']
    )
    
    # Create top layer for ML estimate (red, larger)
    ml_layer = alt.Chart(df_ml).mark_circle(size=100, opacity=1.0).encode(
        x=alt.X('Allele Frequency:Q', title='Allele Frequency'),
        y=alt.Y('Probability:Q', axis=None),
        color=alt.value('red'),
        tooltip=['Allele Frequency', 'Probability', 'Type']
    )
    
    # Layer ML estimate on top
    chart = (base_layer + ml_layer).properties(
        title='Allele Frequency Distribution (ML Estimate in Red)',
        width=500,
        height=300
    ).configure_view(
        strokeWidth=0  # Remove border
    ).configure_axis(
        grid=False  # Remove grid lines
    )
    
    return chart
def visualize_observations(record, sample_name):
    """
    Visualize observations from OBS field
    Two panels: REF allele (left) and ALT allele (right)
    Each panel shows multiple metrics (Posterior Odds, Edit Distance, etc.) as separate stacked bars
    """
    sample = record.samples[sample_name]
    obs = sample['OBS']
    
    # Handle if OBS is a tuple or list
    if isinstance(obs, (tuple, list)):
        obs_string = obs[0] if obs else ""
    else:
        obs_string = obs
    
    # Parse the OBS field
    ref_data = []
    alt_data = []
    
    # Pattern: COUNT + 2-letter-odds + 8 more chars (D,T,A,S,O,P,X,I)
    pattern = r'(\d+)([a-zA-Z]{2})(.{8})'
    matches = re.findall(pattern, obs_string)
    
    for idx, match in enumerate(matches):
        count = int(match[0])
        odds_code = match[1]
        rest = match[2]
        
        # Parse the 8 characters
        edit_distance = rest[0]  # D
        alignment_type = rest[1]  # T
        alt_locus = rest[2]  # A
        strand = rest[3]  # S
        orientation = rest[4]  # O
        read_position = rest[5]  # P
        softclip = rest[6]  # X
        indel = rest[7]  # I
        
        # Determine allele type
        allele_type = odds_code[0].upper()
        kass_raftery = odds_code[1]
        
        # Convert Kass Raftery score to approximate posterior odds
        kr_map = {
            'N': 0.0, 'E': 1.0, 'B': 3.0, 'P': 10.0, 'S': 20.0, 'V': 150.0,
            'n': 0.0, 'e': 0.5, 'b': 1.5, 'p': 5.0, 's': 10.0, 'v': 75.0
        }
        posterior_odds = kr_map.get(kass_raftery, 1.0)
        
        # Get edit distance value
        edit_dist_val = int(edit_distance) if edit_distance.isdigit() else 0
        
        obs_dict = {
            'Observation': str(idx),  # Convert to string for categorical coloring
            'Posterior Odds': posterior_odds,
            'Edit Distance': edit_dist_val,
            'Count': count,
            'Strand': 1 if strand == '+' else (-1 if strand == '-' else 0),
            'Orientation': 1 if orientation == '>' else (-1 if orientation == '<' else 0)
        }
        
        if allele_type == 'A':
            alt_data.append(obs_dict)
        else:
            ref_data.append(obs_dict)
    
    # Create dataframes
    df_ref = pd.DataFrame(ref_data)
    df_alt = pd.DataFrame(alt_data)
    
    # Create REF panel (left)
    if not df_ref.empty:
        # Melt for all metrics
        ref_long = df_ref.melt(
            id_vars=['Observation'],
            value_vars=['Posterior Odds', 'Edit Distance', 'Strand', 'Orientation'],
            var_name='Metric',
            value_name='Value'
        )
        
        ref_chart = alt.Chart(ref_long).mark_bar().encode(
            x=alt.X('Metric:N', title=None, axis=alt.Axis(labelAngle=-45)),
            y=alt.Y('Value:Q', title='Value', stack='zero'),
            color=alt.Color('Observation:N', legend=alt.Legend(title='Observation')),
            order=alt.Order('Observation:N'),
            tooltip=['Observation:N', 'Metric:N', 'Value:Q']
        ).properties(
            title='REF Allele Observations',
            width=400,
            height=400
        )
    else:
        ref_chart = alt.Chart(pd.DataFrame({'text': ['No REF observations']})).mark_text(
            text='No REF observations', size=16
        ).encode().properties(width=400, height=400)
    
    # Create ALT panel (right)
    if not df_alt.empty:
        # Melt for all metrics
        alt_long = df_alt.melt(
            id_vars=['Observation'],
            value_vars=['Posterior Odds', 'Edit Distance', 'Strand', 'Orientation'],
            var_name='Metric',
            value_name='Value'
        )
        
        alt_chart = alt.Chart(alt_long).mark_bar().encode(
            x=alt.X('Metric:N', title=None, axis=alt.Axis(labelAngle=-45)),
            y=alt.Y('Value:Q', title='Value', stack='zero'),
            color=alt.Color('Observation:N', legend=alt.Legend(title='Observation')),
            order=alt.Order('Observation:N'),
            tooltip=['Observation:N', 'Metric:N', 'Value:Q']
        ).properties(
            title='ALT Allele Observations',
            width=400,
            height=400
        )
    else:
        alt_chart = alt.Chart(pd.DataFrame({'text': ['No ALT observations']})).mark_text(
            text='No ALT observations', size=16
        ).encode().properties(width=400, height=400)
    
    # Combine panels side by side
    combined = alt.hconcat(ref_chart, alt_chart)
    
    return combined
# Main execution
if __name__ == "__main__":
    # Open the VCF file
    vcf_file = "examples/example.vcf"
    vcf = pysam.VariantFile(vcf_file)
    
    # Get the first record
    record = next(vcf)
    
    # Get the first sample name
    sample_name = list(record.samples.keys())[0]
    
    print(f"Processing record at {record.chrom}:{record.pos}")
    print(f"Sample: {sample_name}")
    
    # Generate visualizations
    print("\nGenerating event probabilities chart...")
    chart1 = visualize_event_probabilities(record)
    chart1.save('event_probabilities.html')
    print("Saved: event_probabilities.html")
    
    print("\nGenerating allele frequency distribution chart...")
    chart2 = visualize_allele_frequency_distribution(record, sample_name)
    chart2.save('allele_frequency_distribution.html')
    print("Saved: allele_frequency_distribution.html")
    
    print("\nGenerating observations chart...")
    chart3 = visualize_observations(record, sample_name)
    chart3.save('observations.html')
    print("Saved: observations.html")
    
    print("\n✓ All charts generated successfully!")
    
    # Close the VCF file
    vcf.close()

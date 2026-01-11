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
    prob_data = []
    for key, value in record.info.items():
        if key.startswith('PROB_'):
            event_name = key.replace('PROB_', '')
            
            # Handle inf values
            if value == float('inf'):
                probability = 0.0
            else:
                probability = phred_to_prob(value)
            
            prob_data.append({'Event': event_name, 'Probability': probability})
    
    df = pd.DataFrame(prob_data)
    
    print(f"\n=== Event Probabilities ===")
    print(df)
    print(f"Sum of probabilities: {df['Probability'].sum()}")
    
    # Create bar plot WITHOUT log scale
    chart = alt.Chart(df).mark_bar().encode(
        x=alt.X('Event:N', title='Event Type'),
        y=alt.Y('Probability:Q', title='Probability'),  # Removed log scale
        tooltip=['Event', alt.Tooltip('Probability:Q', format='.6f')]
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
    Each metric has a stacked bar with same total height, consistent stacking order
    """
    sample = record.samples[sample_name]
    obs = sample['OBS']
    
    if isinstance(obs, (tuple, list)):
        obs_string = obs[0] if obs else ""
    else:
        obs_string = obs
    
    ref_observations = []
    alt_observations = []
    
    pattern = r'(\d+)([a-zA-Z]{2})(.{8})'
    matches = re.findall(pattern, obs_string)
    
    for idx, match in enumerate(matches):
        count = int(match[0])
        odds_code = match[1]
        rest = match[2]
        
        edit_distance_char = rest[0]
        alignment_type = rest[1]
        alt_locus = rest[2]
        strand = rest[3]
        orientation = rest[4]
        read_position = rest[5]
        softclip = rest[6]
        indel = rest[7]
        
        allele_type = odds_code[0].upper()
        kass_raftery = odds_code[1]
        
        # Map Kass Raftery to descriptive names
        kr_names = {
            'N': 'None', 'E': 'Equal', 'B': 'Barely', 'P': 'Positive', 
            'S': 'Strong', 'V': 'Very Strong',
            'n': 'none', 'e': 'equal', 'b': 'barely', 'p': 'positive',
            's': 'strong', 'v': 'very strong'
        }
        kr_name = kr_names.get(kass_raftery, kass_raftery)
        
        edit_dist_val = int(edit_distance_char) if edit_distance_char.isdigit() else 0
        
        obs_entry = {
            'obs_index': idx,  # Preserve original order for stacking
            'count': count,
            'Posterior Odds': kr_name,  # Category for coloring
            'Edit Distance': str(edit_dist_val),  # Category for coloring
            'Strand': strand,  # +/- for coloring
            'Orientation': orientation,  # >/</*/! for coloring
            'Read Position': read_position,  # ^/* for coloring
            'Softclip': softclip,  # $/. for coloring
            'Indel': indel  # */. for coloring
        }
        
        if allele_type == 'A':
            alt_observations.append(obs_entry)
        else:
            ref_observations.append(obs_entry)
    
    def create_panel(observations, title):
        if not observations:
            return alt.Chart(pd.DataFrame({'text': [f'No {title} observations']})).mark_text(
                text=f'No {title} observations', size=16
            ).encode().properties(width=600, height=400, title=title)
        
        # Create data for each metric
        metrics = ['Posterior Odds', 'Edit Distance', 'Strand', 'Orientation', 
                  'Read Position', 'Softclip', 'Indel']
        
        all_data = []
        for obs in observations:
            obs_index = obs['obs_index']
            count = obs['count']
            
            for metric in metrics:
                category_value = obs[metric]
                
                all_data.append({
                    'obs_index': obs_index,  # For consistent ordering
                    'Metric': metric,
                    'Count': count,
                    'Category': category_value  # The value to color by
                })
        
        df = pd.DataFrame(all_data)
        
        # Create stacked bar chart with different colors per metric
        chart = alt.Chart(df).mark_bar().encode(
            x=alt.X('Metric:N', title=None, axis=alt.Axis(labelAngle=-45)),
            y=alt.Y('Count:Q', title='Count', stack='zero'),
            color=alt.Color('Category:N', 
                          legend=alt.Legend(title='Category'),
                          scale=alt.Scale(scheme='tableau20')),
            order=alt.Order('obs_index:Q'),  # Consistent stacking order
            tooltip=['obs_index:Q', 'Metric:N', 'Category:N', 'Count:Q']
        ).properties(
            title=title,
            width=600,
            height=400
        )
        
        return chart
    
    ref_chart = create_panel(ref_observations, 'REF Allele Observations')
    alt_chart = create_panel(alt_observations, 'ALT Allele Observations')
    
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
    chart1.save('event_probabilities.vl.json')  # Save as Vega-Lite JSON
    print("Saved: event_probabilities.html")
    print("Saved: event_probabilities.vl.json")

    
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

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
    Each metric bar is normalized to 100% showing proportion of each observation
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
        
        kr_map = {
            'N': 0.0, 'E': 1.0, 'B': 3.0, 'P': 10.0, 'S': 20.0, 'V': 150.0,
            'n': 0.0, 'e': 0.5, 'b': 1.5, 'p': 5.0, 's': 10.0, 'v': 75.0
        }
        posterior_odds = kr_map.get(kass_raftery, 1.0)
        
        edit_dist_val = int(edit_distance_char) if edit_distance_char.isdigit() else 0
        
        obs_entry = {
            'obs_id': f'obs_{idx}',
            'count': count,
            'Posterior Odds': posterior_odds,
            'Edit Distance': edit_dist_val,
            'Strand': 1,
            'Orientation': 1,
            'Read Position': 1,
            'Softclip': 1,
            'Indel': 1
        }
        
        if allele_type == 'A':
            alt_observations.append(obs_entry)
        else:
            ref_observations.append(obs_entry)
    
    def create_panel(observations, title):
        if not observations:
            return alt.Chart(pd.DataFrame({'text': [f'No {title} observations']})).mark_text(
                text=f'No {title} observations', size=16
            ).encode().properties(width=500, height=400, title=title)
        
        all_data = []
        for obs in observations:
            obs_id = obs['obs_id']
            count = obs['count']
            
            metrics = ['Posterior Odds', 'Edit Distance', 'Strand', 'Orientation', 
                      'Read Position', 'Softclip', 'Indel']
            
            for metric in metrics:
                if metric in ['Posterior Odds', 'Edit Distance']:
                    value = obs[metric] * count
                else:
                    value = count
                
                all_data.append({
                    'Observation': obs_id,
                    'Metric': metric,
                    'Value': value
                })
        
        df = pd.DataFrame(all_data)
        
        # Use normalize stack to make all bars same height
        chart = alt.Chart(df).mark_bar().encode(
            x=alt.X('Metric:N', title=None, axis=alt.Axis(labelAngle=-45)),
            y=alt.Y('Value:Q', title='Proportion', stack='normalize'),  # Changed to 'normalize'
            color=alt.Color('Observation:N', 
                          legend=alt.Legend(title='Observation'),
                          scale=alt.Scale(scheme='tableau20')),
            order=alt.Order('Observation:N'),
            tooltip=['Observation:N', 'Metric:N', 'Value:Q']
        ).properties(
            title=title,
            width=500,
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

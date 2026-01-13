import altair as alt
import pandas as pd
import pysam
import re

def phred_to_prob(phred_value):
    """Convert PHRED score to probability"""
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
        y=alt.Y('Probability:Q', title='Probability'),
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
    sample = record.samples[sample_name]
    
    afd = sample['AFD']
    af_ml = sample['AF']
    
    if isinstance(af_ml, (tuple, list)):
        af_ml = af_ml[0]
    
    afd_data = []
    
    if isinstance(afd, (tuple, list)):
        afd_entries = afd
    else:
        afd_entries = [afd]
    
    for entry in afd_entries:
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
    
    df_distribution = df[df['Type'] == 'Distribution']
    df_ml = df[df['Type'] == 'ML Estimate']
    
    base_layer = alt.Chart(df_distribution).mark_circle(size=60, opacity=0.7).encode(
        x=alt.X('Allele Frequency:Q', title='Allele Frequency'),
        y=alt.Y('Probability:Q', axis=None),
        color=alt.value('blue'),
        tooltip=['Allele Frequency', 'Probability', 'Type']
    )
    
    ml_layer = alt.Chart(df_ml).mark_circle(size=100, opacity=1.0).encode(
        x=alt.X('Allele Frequency:Q', title='Allele Frequency'),
        y=alt.Y('Probability:Q', axis=None),
        color=alt.value('red'),
        tooltip=['Allele Frequency', 'Probability', 'Type']
    )
    
    chart = (base_layer + ml_layer).properties(
        title='Allele Frequency Distribution (ML Estimate in Red)',
        width=500,
        height=300
    ).configure_view(
        strokeWidth=0
    ).configure_axis(
        grid=False
    )
    
    return chart

def visualize_observations(record, sample_name):
    """
    Visualize observations from OBS field
    Two panels: REF allele (top) and ALT allele (bottom)
    Each metric has its own legend and color meaning
    """
    sample = record.samples[sample_name]
    obs = sample['OBS']

    obs_string = obs[0] if isinstance(obs, (tuple, list)) else obs

    ref_observations = []
    alt_observations = []

    pattern = r'(\d+)([a-zA-Z]{2})(.{8})'
    matches = re.findall(pattern, obs_string)

    for idx, match in enumerate(matches):
        count = int(match[0])
        odds_code = match[1]
        rest = match[2]

        allele_type = odds_code[0].upper()
        kass = odds_code[1]

        kr_map = {
            'N': 'None', 'E': 'Equal', 'B': 'Barely',
            'P': 'Positive', 'S': 'Strong', 'V': 'Very Strong',
            'n': 'none', 'e': 'equal', 'b': 'barely',
            'p': 'positive', 's': 'strong', 'v': 'very strong'
        }

        obs_entry = {
            'Count': count,
            'Posterior Odds': kr_map.get(kass, kass),
            'Edit Distance': rest[0],
            'Strand': rest[3],
            'Orientation': rest[4],
            'Read Position': rest[5],
            'Softclip': rest[6],
            'Indel': rest[7]
        }

        if allele_type == 'A':
            alt_observations.append(obs_entry)
        else:
            ref_observations.append(obs_entry)

    def create_panel(observations, title):
        if not observations:
            return alt.Chart(pd.DataFrame({'x': []})).mark_text(
                text=f'No {title} observations'
            )

        metrics = [
            'Posterior Odds', 'Edit Distance', 'Strand',
            'Orientation', 'Read Position', 'Softclip', 'Indel'
        ]

        charts = []

        for metric in metrics:
            df = pd.DataFrame([
                {'Metric': metric, 'Category': obs[metric], 'Count': obs['Count']}
                for obs in observations
            ])

            chart = alt.Chart(df).mark_bar().encode(
                x=alt.X('Metric:N', axis=alt.Axis(labels=False)),
                y=alt.Y('sum(Count):Q', title='Count'),
                color=alt.Color(
                    'Category:N',
                    title=metric,        # ← LEGEND TITLE IS MEANINGFUL
                    legend=alt.Legend(orient='right')
                ),
                tooltip=['Category:N', 'sum(Count):Q']
            ).properties(
                width=90,
                height=350,
                title=metric
            )

            charts.append(chart)

        return alt.hconcat(*charts).properties(title=title)

    ref_chart = create_panel(ref_observations, 'REF Allele Observations')
    alt_chart = create_panel(alt_observations, 'ALT Allele Observations')

    return alt.vconcat(ref_chart, alt_chart)

# Main execution
if __name__ == "__main__":
    vcf_file = "examples/example.vcf"
    vcf = pysam.VariantFile(vcf_file)
    
    record = next(vcf)
    sample_name = list(record.samples.keys())[0]
    
    print(f"Processing record at {record.chrom}:{record.pos}")
    print(f"Sample: {sample_name}")
    
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
    
    vcf.close()
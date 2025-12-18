"""Core visualization functions for Varlociraptor VCF records."""

import altair as alt
import pandas as pd
from typing import Optional
import pysam


def phred_to_prob(phred_value):
    """Convert PHRED-scaled value to probability."""
    if isinstance(phred_value, (tuple, list)):
        phred_value = phred_value[0]
    return 10 ** (-phred_value / 10)


def visualize_event_probabilities(record: pysam.VariantRecord) -> alt.Chart:
    """Visualize event probabilities from INFO fields starting with PROB_."""
    prob_fields = {}
    for key in record.info:
        if key.startswith('PROB_'):
            event_name = key.replace('PROB_', '')
            phred_value = record.info[key]
            prob_fields[event_name] = phred_to_prob(phred_value)
    
    df = pd.DataFrame([
        {'Event': event, 'Probability': prob}
        for event, prob in prob_fields.items()
    ])
    
    chart = alt.Chart(df).mark_bar().encode(
        x=alt.X('Event:N', title='Event', sort='-y'),
        y=alt.Y('Probability:Q', scale=alt.Scale(type='log'), title='Probability (log scale)'),
        tooltip=['Event', alt.Tooltip('Probability:Q', format='.6f')]
    ).properties(
        title=f'Event Probabilities for {record.chrom}:{record.pos}',
        width=600, height=400
    )
    return chart


def visualize_allele_frequency_distribution(
    record: pysam.VariantRecord, 
    sample_name: Optional[str] = None
) -> alt.Chart:
    """Visualize allele frequency distribution (AFD) for a sample."""
    if sample_name is None:
        sample_name = list(record.samples.keys())[0]
    
    sample = record.samples[sample_name]
    afd_string = sample.get('AFD', '')
    if not afd_string:
        raise ValueError(f"No AFD field found for sample {sample_name}")
    
    afd_data = []
    for entry in afd_string:
        if isinstance(entry, str):
            for pair in entry.split(','):
                if '=' in pair:
                    freq_str, phred_str = pair.split('=')
                    freq = float(freq_str)
                    prob = phred_to_prob(float(phred_str))
                    afd_data.append({'Allele_Frequency': freq, 'Probability': prob})
    
    ml_af = sample.get('AF', None)
    if ml_af is not None:
        if isinstance(ml_af, (list, tuple)):
            ml_af = ml_af[0]
    
    df = pd.DataFrame(afd_data)
    
    if ml_af is not None:
        df['is_ML'] = df['Allele_Frequency'].apply(
            lambda x: 'ML Frequency' if abs(x - ml_af) < 0.001 else 'Other'
        )
    else:
        df['is_ML'] = 'Other'
    
    chart = alt.Chart(df).mark_circle(size=100).encode(
        x=alt.X('Allele_Frequency:Q', title='Allele Frequency', scale=alt.Scale(domain=[0, 1])),
        y=alt.Y('Probability:Q', title='Probability (Density)'),
        color=alt.Color('is_ML:N', scale=alt.Scale(domain=['ML Frequency', 'Other'], range=['red', 'blue'])),
        tooltip=[alt.Tooltip('Allele_Frequency:Q', format='.3f'), alt.Tooltip('Probability:Q', format='.6f'), 'is_ML']
    ).properties(
        title=f'Allele Frequency Distribution - {sample_name} ({record.chrom}:{record.pos})',
        width=600, height=400
    )
    return chart


def visualize_observations(
    record: pysam.VariantRecord, 
    sample_name: Optional[str] = None
) -> alt.Chart:
    """Visualize observations (OBS field)."""
    if sample_name is None:
        sample_name = list(record.samples.keys())[0]
    
    sample = record.samples[sample_name]
    obs_string = sample.get('OBS', '')
    if not obs_string:
        raise ValueError(f"No OBS field found for sample {sample_name}")
    
    df = pd.DataFrame([{'Info': 'OBS data present', 'Value': len(str(obs_string).split(','))}])
    chart = alt.Chart(df).mark_bar().encode(
        x=alt.X('Info:N', title=''),
        y=alt.Y('Value:Q', title='Number of Observations')
    ).properties(
        title=f'Observations - {sample_name} ({record.chrom}:{record.pos})',
        width=400, height=300
    )
    return chart

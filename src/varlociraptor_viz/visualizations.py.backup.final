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
    
    if not prob_fields:
        df = pd.DataFrame({'Event': ['No Data'], 'Probability': [0.5]})
    else:
        df = pd.DataFrame([
            {'Event': event, 'Probability': prob}
            for event, prob in prob_fields.items()
        ])
    
    chart = alt.Chart(df).mark_bar(color='steelblue', size=50).encode(
        x=alt.X('Event:N', title='Event', sort='-y', axis=alt.Axis(labelAngle=-45)),
        y=alt.Y('Probability:Q', scale=alt.Scale(type='log'), title='Probability (log scale)'),
        tooltip=['Event', alt.Tooltip('Probability:Q', format='.6f')]
    ).properties(title=f'Event Probabilities for {record.chrom}:{record.pos}', width=600, height=400)
    
    return chart


def visualize_allele_frequency_distribution(record: pysam.VariantRecord, sample_name: Optional[str] = None) -> alt.Chart:
    """Visualize allele frequency distribution with ML point on top."""
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
    
    if ml_af is not None and not df.empty:
        df['distance'] = df['Allele_Frequency'].apply(lambda x: abs(x - ml_af))
        closest_idx = df['distance'].idxmin()
        df['is_ML'] = 'Other'
        df.at[closest_idx, 'is_ML'] = 'ML'
        df = df.drop('distance', axis=1)
    else:
        df['is_ML'] = 'Other'
    
    df_other = df[df['is_ML'] == 'Other']
    df_ml = df[df['is_ML'] == 'ML']
    
    base = alt.Chart(df_other).mark_circle(size=100, opacity=0.6).encode(
        x=alt.X('Allele_Frequency:Q', title='Allele Frequency', scale=alt.Scale(domain=[0, 1])),
        y=alt.Y('Probability:Q', title=None, axis=None),
        color=alt.value('blue'),
        tooltip=[alt.Tooltip('Allele_Frequency:Q', format='.3f'), alt.Tooltip('Probability:Q', format='.6f')]
    )
    
    if not df_ml.empty:
        ml = alt.Chart(df_ml).mark_circle(size=300, opacity=1.0).encode(
            x='Allele_Frequency:Q', y='Probability:Q', color=alt.value('red'),
            tooltip=[alt.Tooltip('Allele_Frequency:Q', format='.3f', title='ML'), alt.Tooltip('Probability:Q', format='.6f')]
        )
        chart = (base + ml)
    else:
        chart = base
    
    return chart.properties(title=f'AFD - {sample_name} ({record.chrom}:{record.pos})', width=600, height=400).configure_view(strokeWidth=0).configure_axis(grid=False)


def parse_obs_entry(obs_string: str):
    """Parse OBS entry."""
    if not obs_string or len(obs_string) < 10:
        return None
    try:
        parts = list(obs_string)
        return {'count': int(parts[0]), 'direction': parts[1], 'strength': parts[2], 'edit_distance': parts[3],
                'alignment_type': parts[4], 'alt_locus': parts[5], 'strand': parts[6], 'orientation': parts[7],
                'read_position': parts[8], 'softclip': parts[9], 'indel': parts[10] if len(parts) > 10 else '.'}
    except (ValueError, IndexError):
        return None


def get_posterior_odds_value(direction: str, strength: str) -> float:
    """Convert Kass-Raftery score."""
    strength_upper = strength.upper()
    odds_map = {'N': 0.1, 'E': 1, 'B': 5, 'P': 30, 'S': 300, 'V': 3000}
    base_odds = odds_map.get(strength_upper, 1)
    if direction == 'R':
        base_odds = -base_odds
    if strength.islower():
        base_odds *= 0.7
    return base_odds


def visualize_observations(record: pysam.VariantRecord, sample_name: Optional[str] = None) -> alt.Chart:
    """Visualize observations with REF/ALT separation."""
    if sample_name is None:
        sample_name = list(record.samples.keys())[0]
    
    sample = record.samples[sample_name]
    obs_string = sample.get('OBS', '')
    if not obs_string:
        raise ValueError(f"No OBS field found for sample {sample_name}")
    
    obs_entries = obs_string if isinstance(obs_string, (list, tuple)) else str(obs_string).split(',')
    
    ref_data, alt_data = [], []
    for entry in obs_entries:
        parsed = parse_obs_entry(entry)
        if parsed:
            (alt_data if parsed['direction'] == 'A' else ref_data).append(parsed)
    
    def create_panel(data, allele):
        rows = []
        for i, obs in enumerate(data):
            obs_id = f"{allele}_{i}"
            rows.append({'Observation': obs_id, 'Attribute': 'Posterior Odds', 'Value': abs(get_posterior_odds_value(obs['direction'], obs['strength'])),
                        'Count': obs['count'], 'Category': f"{obs['direction']}{obs['strength']}", 'Allele': allele})
            if obs['edit_distance'] != '.':
                try:
                    rows.append({'Observation': obs_id, 'Attribute': 'Edit Distance', 'Value': int(obs['edit_distance']),
                                'Count': obs['count'], 'Category': obs['edit_distance'], 'Allele': allele})
                except ValueError:
                    pass
            rows.append({'Observation': obs_id, 'Attribute': 'Strand', 'Value': {'+': 1, '-': 1, '*': 1, '.': 0}.get(obs['strand'], 0),
                        'Count': obs['count'], 'Category': obs['strand'], 'Allele': allele})
            rows.append({'Observation': obs_id, 'Attribute': 'Softclip', 'Value': 1 if obs['softclip'] == '$' else 0,
                        'Count': obs['count'], 'Category': obs['softclip'], 'Allele': allele})
            rows.append({'Observation': obs_id, 'Attribute': 'Indel', 'Value': 1 if obs['indel'] == '*' else 0,
                        'Count': obs['count'], 'Category': obs['indel'], 'Allele': allele})
        return pd.DataFrame(rows)
    
    df = pd.concat([create_panel(ref_data, 'REF'), create_panel(alt_data, 'ALT')], ignore_index=True)
    if df.empty:
        return alt.Chart(pd.DataFrame()).mark_bar()
    
    return alt.Chart(df).mark_bar().encode(
        x=alt.X('Observation:N', title='Observation', axis=alt.Axis(labels=False)),
        y=alt.Y('Value:Q', title='Value'), color=alt.Color('Category:N', title='Category'),
        column=alt.Column('Attribute:N', title=None), row=alt.Row('Allele:N', title='Allele Type'),
        tooltip=['Observation', 'Attribute', 'Value', 'Count', 'Category', 'Allele']
    ).properties(width=150, height=200, title=f'Observations - {sample_name} ({record.chrom}:{record.pos})').resolve_scale(y='independent')

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
        df = pd.DataFrame({'Event': ['No PROB fields'], 'Probability': [0.5]})
    else:
        df = pd.DataFrame([
            {'Event': event, 'Probability': prob}
            for event, prob in prob_fields.items()
        ])
    
    chart = alt.Chart(df).mark_bar(
        color='steelblue',
        size=50
    ).encode(
        x=alt.X('Event:N', title='Event', sort='-y', axis=alt.Axis(labelAngle=-45)),
        y=alt.Y('Probability:Q', 
                scale=alt.Scale(type='log'),
                title='Probability (log scale)'),
        tooltip=['Event', alt.Tooltip('Probability:Q', format='.6f')]
    ).properties(
        title=f'Event Probabilities for {record.chrom}:{record.pos}',
        width=600,
        height=400
    )
    
    return chart


def visualize_allele_frequency_distribution(
    record: pysam.VariantRecord, 
    sample_name: Optional[str] = None
) -> alt.Chart:
    """
    Visualize allele frequency distribution (AFD) for a sample.
    
    Professor's requirements:
    1. ML point displayed ON TOP using layering
    2. Y-axis removed (density values not meaningful)
    3. Grid removed
    """
    if sample_name is None:
        sample_name = list(record.samples.keys())[0]
    
    sample = record.samples[sample_name]
    afd_string = sample.get('AFD', '')
    if not afd_string:
        raise ValueError(f"No AFD field found for sample {sample_name}")
    
    # Parse AFD field
    afd_data = []
    for entry in afd_string:
        if isinstance(entry, str):
            for pair in entry.split(','):
                if '=' in pair:
                    freq_str, phred_str = pair.split('=')
                    freq = float(freq_str)
                    prob = phred_to_prob(float(phred_str))
                    afd_data.append({
                        'Allele_Frequency': freq, 
                        'Probability': prob
                    })
    
    # Get ML frequency
    ml_af = sample.get('AF', None)
    if ml_af is not None:
        if isinstance(ml_af, (list, tuple)):
            ml_af = ml_af[0]
    
    df = pd.DataFrame(afd_data)
    
    # Identify ML point (find closest to AF value)
    if ml_af is not None and not df.empty:
        df['is_ML'] = df['Allele_Frequency'].apply(
            lambda x: abs(x - ml_af) < 0.01  # Tolerance of 0.01
        )
    else:
        df['is_ML'] = False
    
    # Separate into two dataframes for layering
    df_other = df[df['is_ML'] == False].copy()
    df_ml = df[df['is_ML'] == True].copy()
    
    # LAYER 1: Base layer with non-ML points (blue, smaller, drawn first)
    base_layer = alt.Chart(df_other).mark_circle(
        size=100,
        color='blue',
        opacity=0.6,
        filled=True
    ).encode(
        x=alt.X('Allele_Frequency:Q', 
                title='Allele Frequency',
                scale=alt.Scale(domain=[0, 1]),
                axis=alt.Axis(grid=False)),
        y=alt.Y('Probability:Q', 
                axis=None),  # NO Y-AXIS as requested
        tooltip=[
            alt.Tooltip('Allele_Frequency:Q', format='.3f', title='Frequency'),
            alt.Tooltip('Probability:Q', format='.6f', title='Probability')
        ]
    )
    
    # LAYER 2: ML layer (red, larger, drawn ON TOP)
    if not df_ml.empty:
        ml_layer = alt.Chart(df_ml).mark_circle(
            size=400,      # Much larger
            color='red',
            opacity=1.0,
            filled=True,
            stroke='darkred',
            strokeWidth=2
        ).encode(
            x=alt.X('Allele_Frequency:Q'),
            y=alt.Y('Probability:Q'),
            tooltip=[
                alt.Tooltip('Allele_Frequency:Q', format='.3f', title='ML Frequency'),
                alt.Tooltip('Probability:Q', format='.6f', title='Probability')
            ]
        )
        
        # SUPERPOSITION: Layer ML on top of base
        combined_chart = (base_layer + ml_layer)
    else:
        combined_chart = base_layer
    
    # Configure chart (remove grid as requested)
    final_chart = combined_chart.properties(
        title=f'Allele Frequency Distribution - {sample_name} ({record.chrom}:{record.pos})',
        width=600,
        height=400
    ).configure_view(
        strokeWidth=0  # No border
    ).configure_axis(
        grid=False     # NO GRID as requested
    )
    
    return final_chart


def parse_obs_entry(obs_string: str):
    """
    Parse a single OBS field entry.
    Format: CBDTASOPXI
    """
    if not obs_string or len(obs_string) < 10:
        return None
    
    try:
        parts = list(obs_string)
        return {
            'count': int(parts[0]),
            'direction': parts[1],      # 'A' = ALT, 'R' = REF
            'strength': parts[2],        # B, P, S, V (Kass-Raftery)
            'edit_distance': parts[3],
            'alignment_type': parts[4],
            'alt_locus': parts[5],
            'strand': parts[6],          # +, -, *
            'orientation': parts[7],
            'read_position': parts[8],
            'softclip': parts[9],        # $, .
            'indel': parts[10] if len(parts) > 10 else '.'
        }
    except (ValueError, IndexError):
        return None


def get_posterior_odds_value(direction: str, strength: str) -> float:
    """Convert Kass-Raftery score to numeric value."""
    strength_upper = strength.upper()
    
    # Kass-Raftery mapping
    odds_map = {
        'N': 0.1,    # None
        'E': 1,      # Equal
        'B': 5,      # Barely
        'P': 30,     # Positive
        'S': 300,    # Strong
        'V': 3000    # Very strong
    }
    
    base_odds = odds_map.get(strength_upper, 1)
    
    # Direction affects sign
    if direction == 'R':
        base_odds = -base_odds
    
    # Lowercase indicates low mapping quality
    if strength.islower():
        base_odds *= 0.7
    
    return base_odds


def visualize_observations(
    record: pysam.VariantRecord, 
    sample_name: Optional[str] = None
) -> alt.Chart:
    """
    Visualize observations (OBS field) separated by REF and ALT alleles.
    
    Professor's requirements:
    1. Multiple bars (not just count)
    2. Coloring by category
    3. Distinguish between REF and ALT
    """
    if sample_name is None:
        sample_name = list(record.samples.keys())[0]
    
    sample = record.samples[sample_name]
    obs_string = sample.get('OBS', '')
    if not obs_string:
        raise ValueError(f"No OBS field found for sample {sample_name}")
    
    # Split observations
    if isinstance(obs_string, (list, tuple)):
        obs_entries = obs_string
    else:
        obs_entries = str(obs_string).split(',')
    
    # Separate REF and ALT observations
    ref_observations = []
    alt_observations = []
    
    for entry in obs_entries:
        parsed = parse_obs_entry(entry)
        if parsed:
            if parsed['direction'] == 'A':
                alt_observations.append(parsed)
            else:
                ref_observations.append(parsed)
    
    # Create data for visualization
    def create_viz_data(observations, allele_type):
        """Create rows for multiple bar attributes."""
        rows = []
        
        for i, obs in enumerate(observations):
            obs_id = f"{allele_type}_{i+1}"
            
            # Attribute 1: Posterior Odds
            odds = get_posterior_odds_value(obs['direction'], obs['strength'])
            rows.append({
                'Observation_ID': obs_id,
                'Allele': allele_type,
                'Attribute': 'Posterior_Odds',
                'Value': abs(odds),
                'Category': f"{obs['direction']}{obs['strength']}"
            })
            
            # Attribute 2: Edit Distance
            if obs['edit_distance'] not in ['.', '*']:
                try:
                    edit_dist = int(obs['edit_distance'])
                    rows.append({
                        'Observation_ID': obs_id,
                        'Allele': allele_type,
                        'Attribute': 'Edit_Distance',
                        'Value': edit_dist,
                        'Category': f"ED_{edit_dist}"
                    })
                except ValueError:
                    pass
            
            # Attribute 3: Strand
            strand_value = 1 if obs['strand'] in ['+', '-'] else 0
            rows.append({
                'Observation_ID': obs_id,
                'Allele': allele_type,
                'Attribute': 'Strand',
                'Value': strand_value,
                'Category': f"Strand_{obs['strand']}"
            })
            
            # Attribute 4: Softclip
            softclip_value = 1 if obs['softclip'] == '$' else 0
            rows.append({
                'Observation_ID': obs_id,
                'Allele': allele_type,
                'Attribute': 'Softclip',
                'Value': softclip_value,
                'Category': f"Softclip_{obs['softclip']}"
            })
            
            # Attribute 5: Indel
            indel_value = 1 if obs['indel'] == '*' else 0
            rows.append({
                'Observation_ID': obs_id,
                'Allele': allele_type,
                'Attribute': 'Indel',
                'Value': indel_value,
                'Category': f"Indel_{obs['indel']}"
            })
        
        return rows
    
    # Combine REF and ALT data
    all_rows = []
    all_rows.extend(create_viz_data(ref_observations, 'REF'))
    all_rows.extend(create_viz_data(alt_observations, 'ALT'))
    
    df = pd.DataFrame(all_rows)
    
    if df.empty:
        return alt.Chart(pd.DataFrame()).mark_bar()
    
    # Create faceted chart with COLORING and REF/ALT separation
    chart = alt.Chart(df).mark_bar().encode(
        x=alt.X('Observation_ID:N', 
                title='Observation',
                axis=alt.Axis(labels=False, ticks=False)),
        y=alt.Y('Value:Q', 
                title='Value'),
        color=alt.Color('Category:N',           # COLORING by category
                       title='Category',
                       legend=alt.Legend(orient='right')),
        column=alt.Column('Attribute:N',        # MULTIPLE BARS (columns)
                         title='Attribute Type'),
        row=alt.Row('Allele:N',                 # REF vs ALT (rows)
                   title='Allele Type'),
        tooltip=['Observation_ID', 'Allele', 'Attribute', 'Value', 'Category']
    ).properties(
        width=120,
        height=180,
        title=f'Observations - {sample_name} ({record.chrom}:{record.pos})'
    ).resolve_scale(
        y='independent'  # Each attribute has own scale
    )
    
    return chart

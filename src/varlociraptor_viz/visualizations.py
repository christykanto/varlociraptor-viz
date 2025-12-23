import altair as alt
import pandas as pd
from typing import Optional
import pysam


def phred_to_prob(phred_value):
    if isinstance(phred_value, (tuple, list)):
        phred_value = phred_value[0]
    return 10 ** (-phred_value / 10)


def visualize_event_probabilities(record: pysam.VariantRecord) -> alt.Chart:
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
id removed
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
    
    df_other = df[df['is_ML'] == False].copy()
    df_ml = df[df['is_ML'] == True].copy()
    
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
                axis=None),  
        tooltip=[
            alt.Tooltip('Allele_Frequency:Q', format='.3f', title='Frequency'),
            alt.Tooltip('Probability:Q', format='.6f', title='Probability')
        ]
    )
    
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
        
        combined_chart = (base_layer + ml_layer)
    else:
        combined_chart = base_layer
    
    final_chart = combined_chart.properties(
        title=f'Allele Frequency Distribution - {sample_name} ({record.chrom}:{record.pos})',
        width=600,
        height=400
    ).configure_view(
        strokeWidth=0 
    ).configure_axis(
        grid=False     
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
            'direction': parts[1],     
            'strength': parts[2],      
            'edit_distance': parts[3],
            'alignment_type': parts[4],
            'alt_locus': parts[5],
            'strand': parts[6],         
            'orientation': parts[7],
            'read_position': parts[8],
            'softclip': parts[9],       
            'indel': parts[10] if len(parts) > 10 else '.'
        }
    except (ValueError, IndexError):
        return None


def get_posterior_odds_value(direction: str, strength: str) -> float:
    strength_upper = strength.upper()
    
    odds_map = {
        'N': 0.1,    
        'E': 1,      
        'B': 5,      
        'P': 30,     
        'S': 300,    
        'V': 3000    
    }
    
    base_odds = odds_map.get(strength_upper, 1)
    
   
    if direction == 'R':
        base_odds = -base_odds
    
    if strength.islower():
        base_odds *= 0.7
    
    return base_odds


def visualize_observations(
    record: pysam.VariantRecord, 
    sample_name: Optional[str] = None
) -> alt.Chart:
   
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
        rows = []
        
        for i, obs in enumerate(observations):
            obs_id = f"{allele_type}_{i+1}"
            
            odds = get_posterior_odds_value(obs['direction'], obs['strength'])
            rows.append({
                'Observation_ID': obs_id,
                'Allele': allele_type,
                'Attribute': 'Posterior_Odds',
                'Value': abs(odds),
                'Category': f"{obs['direction']}{obs['strength']}"
            })
            
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
            
            strand_value = 1 if obs['strand'] in ['+', '-'] else 0
            rows.append({
                'Observation_ID': obs_id,
                'Allele': allele_type,
                'Attribute': 'Strand',
                'Value': strand_value,
                'Category': f"Strand_{obs['strand']}"
            })
            
            softclip_value = 1 if obs['softclip'] == '$' else 0
            rows.append({
                'Observation_ID': obs_id,
                'Allele': allele_type,
                'Attribute': 'Softclip',
                'Value': softclip_value,
                'Category': f"Softclip_{obs['softclip']}"
            })
            
            indel_value = 1 if obs['indel'] == '*' else 0
            rows.append({
                'Observation_ID': obs_id,
                'Allele': allele_type,
                'Attribute': 'Indel',
                'Value': indel_value,
                'Category': f"Indel_{obs['indel']}"
            })
        
        return rows
    
    all_rows = []
    all_rows.extend(create_viz_data(ref_observations, 'REF'))
    all_rows.extend(create_viz_data(alt_observations, 'ALT'))
    
    df = pd.DataFrame(all_rows)
    
    if df.empty:
        return alt.Chart(pd.DataFrame()).mark_bar()
    
    chart = alt.Chart(df).mark_bar().encode(
        x=alt.X('Observation_ID:N', 
                title='Observation',
                axis=alt.Axis(labels=False, ticks=False)),
        y=alt.Y('Value:Q', 
                title='Value'),
        color=alt.Color('Category:N',           
                       title='Category',
                       legend=alt.Legend(orient='right')),
        column=alt.Column('Attribute:N',       
                         title='Attribute Type'),
        row=alt.Row('Allele:N',                
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

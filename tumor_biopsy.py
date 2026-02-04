import pysam
from visualize_vcf import visualize_observations, visualize_allele_frequency_distribution

bcf_file = "examples/calls.bcf"
vcf = pysam.VariantFile(bcf_file)

record = next(vcf)

sample_name = "tumor_biopsy"

print(f"Processing record at {record.chrom}:{record.pos}")
print(f"Sample: {sample_name}")
print(f"Available samples: {list(record.samples.keys())}")

print("\nGenerating observations chart for tumor_biopsy...")
chart_obs = visualize_observations(record, sample_name)
chart_obs.save('observations_tumor_biopsy.html')
print("Saved: observations_tumor_biopsy.html")

print("\nGenerating allele frequency distribution chart for tumor_biopsy...")
chart_afd = visualize_allele_frequency_distribution(record, sample_name)
chart_afd.save('afd_tumor_biopsy.html')
print("Saved: afd_tumor_biopsy.html")

print("\n✓ Charts generated successfully!")

vcf.close()
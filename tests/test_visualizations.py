import pytest
import pysam
from varlociraptor_viz import (
    visualize_event_probabilities,
    visualize_allele_frequency_distribution,
    visualize_observations,
    phred_to_prob
)


def test_phred_to_prob():
    assert abs(phred_to_prob(0) - 1.0) < 0.001
    assert abs(phred_to_prob(10) - 0.1) < 0.001
    assert abs(phred_to_prob(20) - 0.01) < 0.001


def test_phred_to_prob_tuple():

    assert abs(phred_to_prob((10, 20)) - 0.1) < 0.001


@pytest.fixture
def example_vcf_path(tmp_path):
    vcf_content = """##fileformat=VCFv4.2
##INFO=<ID=PROB_SOMATIC,Number=1,Type=Float,Description="Somatic probability">
##INFO=<ID=PROB_ABSENT,Number=1,Type=Float,Description="Absent probability">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AF,Number=1,Type=Float,Description="Allele frequency">
##FORMAT=<ID=AFD,Number=.,Type=String,Description="Allele frequency distribution">
##FORMAT=<ID=OBS,Number=.,Type=String,Description="Observations">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1
chr1\t1000\t.\tA\tT\t.\t.\tPROB_SOMATIC=5.2;PROB_ABSENT=30.0\tGT:AF:AFD:DP:OBS\t0/1:0.42:0.42=0.5,0.5=3.2:50:1AB3sp+>.$.,2AS0p.+>^.$*
"""
    vcf_file = tmp_path / "test.vcf"
    vcf_file.write_text(vcf_content)
    return str(vcf_file)


def test_visualize_event_probabilities(example_vcf_path):
    vcf = pysam.VariantFile(example_vcf_path)
    record = next(vcf)
    chart = visualize_event_probabilities(record)
    assert chart is not None
    assert chart.mark == 'bar'
    vcf.close()


def test_visualize_allele_frequency_distribution(example_vcf_path):
    vcf = pysam.VariantFile(example_vcf_path)
    record = next(vcf)
    chart = visualize_allele_frequency_distribution(record, "sample1")
    assert chart is not None
    vcf.close()


def test_visualize_observations(example_vcf_path):
    vcf = pysam.VariantFile(example_vcf_path)
    record = next(vcf)
    chart = visualize_observations(record, "sample1")
    assert chart is not None
    vcf.close()

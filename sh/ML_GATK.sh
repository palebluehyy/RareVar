REF=$1
BAM=$2
VAR=$3
BED=$4
NT=$5

set -e

echo_step()
{
	str=$1
	echo "####################"
	echo "# $str"
	echo "####################"
}

if [ ! -s ${VAR%.in.vcf}.vcf ]; then
	echo_step "Adding sequencing quality annotation..."
	java -Xmx4G -jar /home/bin/GenomeAnalysisTK.jar \
	-nt $NT \
	-R $REF \
	-T VariantAnnotator \
	-I $BAM \
	-o ${VAR%.in.vcf}.vcf \
	-dt NONE \
	-A GCContent \
	-A HomopolymerRun \
	-A SpanningDeletions \
	-A BaseQualityRankSumTest \
	-A FisherStrand \
	-A HaplotypeScore \
	-A MappingQualityRankSumTest \
	-A QualByDepth \
	-A ReadPosRankSumTest \
	-A RMSMappingQuality \
	-A StrandOddsRatio \
	-L $BED \
	--variant $VAR \
	--allow_potentially_misencoded_quality_scores \
	-U ALL
fi

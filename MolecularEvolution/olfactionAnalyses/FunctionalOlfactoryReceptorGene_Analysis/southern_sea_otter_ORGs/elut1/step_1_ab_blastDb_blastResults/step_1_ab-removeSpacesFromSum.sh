# remove spaces from .sum output:
spp=elut2
sed 's/\/ /\//g' olfac_5sp_2seq_query.$spp.sum > olfac_5sp_2seq_query.$spp.NOSPACES.sum
# keep the space in the query name (e.g. Query_1; MamuOR14.5.7_CladeA not Query_1;MamuOR14.5.7_CladeA)
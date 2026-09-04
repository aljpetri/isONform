# Test data

`sirv_sim/` is a small dataset for checking that an installation works. Two
clusters, 51 reads each, in the shape `isONform_parallel --fastq_folder` expects:
one fastq per cluster, named by cluster id.

## Provenance

Simulated SIRV reads (7% error, the `reads_10k_err7%` set), passed through the
same upstream pipeline every corpus here uses --- isONclust, then **isONcorrect**
--- and then subsampled to the three most abundant transcripts of two clusters,
17 reads each.

**The reads are corrected, and that matters.** isONform consumes corrected reads:
on raw simulated reads at 7% error no two reads take the same route through the
graph, every read becomes its own isoform, and the run tells you nothing. A first
version of this test set was built from raw reads and produced 102 isoforms from
102 reads --- which is what that failure looks like, and why it is written down.

| cluster | reads | transcripts |
|---|---|---|
| `0` | 51 | SIRV101, SIRV102, SIRV103 |
| `1` | 51 | SIRV201, SIRV202, SIRV204 |

## Expected result

```
isONform_parallel --fastq_folder test_data/sirv_sim --outfolder /tmp/isonform_test \
                  --t 4 --split_wrt_batches --iso_abundance 3
```

**7 isoforms**, recovering all six transcripts above at >= 0.993 identity, with
SIRV201 emitted twice. The counts are stable: `--iso_abundance 1` gives 12 and
`5` gives 6.

Exact sequences are not promised, and `--faithful` may differ from the default on
a set this small --- a handful of reads is enough for one alignment verdict to
change an isoform. Judge the run by the transcripts recovered, not by a checksum.

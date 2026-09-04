//! Reading fastq and fasta, matching `help_functions.readfq`.
//!
//! The reference uses Heng Li's `readfq` generator. Three of its behaviours reach
//! the output and are reproduced here:
//!
//! * **spaces in the header become underscores** (`last[1:].replace(" ", "_")`),
//!   and the name is the *whole* header line, not the part before the first
//!   space. Read accessions appear verbatim in `mapping.txt` and the skip file,
//!   so this is observable;
//! * a record whose quality line is shorter than its sequence is yielded as
//!   **fasta** (quality `None`) rather than raising;
//! * `@` inside a quality string is not a record boundary --- the quality loop
//!   reads by *length*, not by looking for the next header, which is what makes
//!   fastq parseable at all.

/// One record: accession, sequence, and quality when the file had one.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Record {
    pub acc: String,
    pub seq: Vec<u8>,
    pub qual: Option<Vec<u8>>,
}

/// `readfq`, over the lines of a file.
///
/// Line endings are stripped as the reference's `l[:-1]` does. A trailing `\r`
/// from a CRLF file is *not* stripped by the reference, so it is not stripped
/// here either --- it would end up inside the sequence there too.
pub fn read_fastq(text: &str) -> Vec<Record> {
    let mut out = Vec::new();
    // `l[:-1]` drops the last character of every line, which for a file whose
    // final line has no newline drops a real character. Splitting on '\n' and
    // keeping what is left reproduces the common case and is kinder to the
    // pathological one.
    let mut lines = text.split('\n').peekable();
    let mut last: Option<String> = None;

    loop {
        if last.is_none() {
            // Search for the next header.
            for l in lines.by_ref() {
                if l.starts_with('>') || l.starts_with('@') {
                    last = Some(l.to_string());
                    break;
                }
            }
        }
        let Some(header) = last.take() else { break };
        let acc = header[1..].replace(' ', "_");

        let mut seqs: Vec<&str> = Vec::new();
        let mut boundary: Option<String> = None;
        for l in lines.by_ref() {
            if l.starts_with('@') || l.starts_with('+') || l.starts_with('>') {
                boundary = Some(l.to_string());
                break;
            }
            seqs.push(l);
        }
        let seq: String = seqs.concat();

        match &boundary {
            // Not a `+` line: this was a fasta record.
            None | Some(_) if boundary.as_ref().is_none_or(|b| !b.starts_with('+')) => {
                out.push(Record {
                    acc,
                    seq: seq.into_bytes(),
                    qual: None,
                });
                last = boundary;
                if last.is_none() {
                    break;
                }
            }
            _ => {
                // fastq: read quality until it is at least as long as the seq.
                let mut quals: Vec<&str> = Vec::new();
                let mut len = 0usize;
                let mut complete = false;
                for l in lines.by_ref() {
                    quals.push(l);
                    len += l.len();
                    if len >= seq.len() {
                        complete = true;
                        break;
                    }
                }
                if complete {
                    out.push(Record {
                        acc,
                        seq: seq.into_bytes(),
                        qual: Some(quals.concat().into_bytes()),
                    });
                    last = None;
                } else {
                    // EOF before enough quality: yielded as fasta, then stop.
                    out.push(Record {
                        acc,
                        seq: seq.into_bytes(),
                        qual: None,
                    });
                    break;
                }
            }
        }
        if lines.peek().is_none() && last.is_none() {
            break;
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn accs(rs: &[Record]) -> Vec<&str> {
        rs.iter().map(|r| r.acc.as_str()).collect()
    }

    #[test]
    fn reads_a_plain_fastq() {
        let got = read_fastq("@r1\nACGT\n+\nIIII\n@r2\nTTTT\n+\nJJJJ\n");
        assert_eq!(accs(&got), vec!["r1", "r2"]);
        assert_eq!(got[0].seq, b"ACGT");
        assert_eq!(got[0].qual.as_deref(), Some(&b"IIII"[..]));
        assert_eq!(got[1].seq, b"TTTT");
    }

    #[test]
    fn spaces_in_the_header_become_underscores() {
        // `last[1:].replace(" ", "_")`, and the *whole* header is kept --- not
        // truncated at the first space, which is what most parsers do. The
        // accession reaches `mapping.txt` verbatim, so this is observable.
        let got = read_fastq("@read 1 extra info\nACGT\n+\nIIII\n");
        assert_eq!(got[0].acc, "read_1_extra_info");
    }

    #[test]
    fn an_at_sign_inside_the_quality_is_not_a_record_boundary() {
        // The quality loop reads by length, so a quality string starting with
        // '@' does not start a new record. Splitting on headers instead would
        // silently lose reads.
        let got = read_fastq("@r1\nACGTACGT\n+\n@@@@@@@@\n@r2\nTTTT\n+\nIIII\n");
        assert_eq!(accs(&got), vec!["r1", "r2"]);
        assert_eq!(got[0].qual.as_deref(), Some(&b"@@@@@@@@"[..]));
    }

    #[test]
    fn a_fasta_record_has_no_quality() {
        let got = read_fastq(">r1\nACGT\n>r2\nTTTT\n");
        assert_eq!(accs(&got), vec!["r1", "r2"]);
        assert!(got.iter().all(|r| r.qual.is_none()));
    }

    #[test]
    fn a_multiline_sequence_is_joined() {
        let got = read_fastq(">r1\nACGT\nACGT\nACGT\n");
        assert_eq!(got[0].seq, b"ACGTACGTACGT");
    }

    #[test]
    fn a_truncated_quality_yields_a_fasta_record_rather_than_raising() {
        // The reference's `if last:` fallback after the quality loop.
        let got = read_fastq("@r1\nACGTACGT\n+\nII\n");
        assert_eq!(got.len(), 1);
        assert_eq!(got[0].seq, b"ACGTACGT");
        assert_eq!(got[0].qual, None, "yielded as fasta");
    }

    #[test]
    fn an_empty_file_yields_nothing() {
        assert!(read_fastq("").is_empty());
        assert!(read_fastq("\n\n").is_empty());
    }
}

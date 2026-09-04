# Window-length splice test + report figure generation

[← Back to repo README](../../README.md) | See also: [docs/round2_hypothesis_evaluation.tex](../../docs/round2_hypothesis_evaluation.pdf) | [docs/notebook5_revamp_progress.md](../../docs/notebook5_revamp_progress.md) (2026-09-04 log)

Two scripts, run on bluehive (raw `.mat` data lives there), supporting
`docs/round2_hypothesis_evaluation.tex`:

## `splice_test.py`

Tests the window-length recommendation from the report's Section 4.2 directly: splices several
genuinely non-overlapping, contiguous 3-hour windows (confirmed from
`python/ccf_pipeline/prepare_data.py`'s 50%-overlap windowing formula -- windows two apart are
exactly contiguous, no time duplicated or skipped) into a longer continuous trace, and compares
FastMspec+picker convergence at the standard 3-hour window vs. the spliced ~24-hour one, for
`AFSKRH_XVBAEL` (quartile 4's mean-distance pair). **Result: did not converge either way**,
even in the corrected version that preserves the full `coh_num=107` day stack (an earlier attempt
using only 10 days confounded the test by collapsing `coh_num` to 10). Reported as a genuine
negative result in the report, not smoothed over -- see its Section 4.3.5. `N_DAYS_TO_SPLICE`
controls how many days to splice (`None` = all available).

## `gen_report_figures.py`

Generates the four coherence-spectrum + KDE/picker-diagnostic figures used in the report's
Section 5 (one representative pair per distance quartile, real data, `plotting=True` on the
vendored picker). Edit `EXAMPLES` to regenerate for different pairs/bandwidths.

## Running

```bash
scp splice_test.py gen_report_figures.py bluehive:/scratch/tolugboj_lab/FastMSPEC_dispcurve_batch/
ssh bluehive
source /scratch/tolugboj_lab/softwares/anaconda/anaconda3/2021.05/etc/profile.d/conda.sh
conda activate fastmspec_batch
cd /scratch/tolugboj_lab/FastMSPEC_dispcurve_batch
python3 splice_test.py          # ~30-90 min, submit via sbatch on a real run
python3 gen_report_figures.py   # ~1-1.5hr for all 4 examples, submit via sbatch
```
Raw `.mat` data (not committed, matches this repo's established practice) is read via the
standard manifest (`data/madagascar_stn_conn_ccflist.csv`'s `filelocation` column).

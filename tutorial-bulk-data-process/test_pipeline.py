import unittest
import numpy as np
import pandas as pd
from pipeline import rank_score, choose_tcga, build_os, stage_number


class PipelineTests(unittest.TestCase):
    def test_nominal_threshold_does_not_use_fdr_or_effect_size_cutoff(self):
        from canellas_threshold_analysis import nominal_pm_up
        results = pd.DataFrame({"pvalue": [.049,.05,.001,.001,.01,np.nan],
            "log2FoldChange": [.01,2,0,-1,.6,1], "padj": [.5,.05,.01,.01,.8,np.nan]},
            index=["small_effect", "boundary", "zero", "negative", "large_fdr", "missing"])
        self.assertEqual(results.index[nominal_pm_up(results)].tolist(), ["small_effect", "large_fdr"])

    def test_rank_score_and_order_invariance(self):
        expression = pd.DataFrame({"a": [0, 1, 1, 3], "b": [3, 2, 1, 0]}, index=list("ABCD"))
        score, coverage = rank_score(expression, ["B", "D"], list("ABCD"))
        np.testing.assert_allclose(score, [.8125, .5])
        shuffled, _ = rank_score(expression.loc[list("DCBA"), ["b", "a"]], ["D", "B"], list("ABCD"))
        pd.testing.assert_series_equal(score.sort_index(), shuffled.sort_index())
        self.assertEqual(coverage, 1)

    def test_coverage_and_empty_signature(self):
        x = pd.DataFrame({"sample": [1, 2, 3, 4]}, index=list("ABCD"))
        _, coverage = rank_score(x, list("ABCDE"), list("ABCD"))
        self.assertEqual(coverage, .8)
        with self.assertRaises(ValueError):
            rank_score(x, list("ABCDEF"), list("ABCD"))
        with self.assertRaises(ValueError):
            rank_score(x, [], list("ABCD"))

    def test_tcga_selection(self):
        d = pd.DataFrame({"Case ID": ["TCGA-AA-0001"]*4,
            "Sample ID": ["TCGA-AA-0001-01B", "TCGA-AA-0001-01A", "TCGA-AA-0001-01A", "TCGA-AA-0001-11A"],
            "File ID": ["a", "z", "b", "c"], "Tumor Descriptor": ["Primary"]*3+["Not Applicable"]})
        audit, selected = choose_tcga(d)
        self.assertEqual(selected["File ID"].tolist(), ["b"])
        self.assertEqual(audit.selection_reason.eq("selected").sum(), 1)

    def test_os_dedup_and_invalid_death(self):
        d = pd.DataFrame({"cases.submitter_id": ["p1", "p1", "p2", "p3", "p4"],
            "diagnoses.diagnosis_is_primary_disease": [True]*5, "diagnoses.days_to_diagnosis": [0]*5,
            "demographic.vital_status": ["Alive", "Alive", "Dead", "Dead", "Dead"],
            "demographic.days_to_death": ["'--", "'--", "'--", 400, 0],
            "diagnoses.days_to_last_follow_up": [100,100,200,300,0],
            "diagnoses.age_at_diagnosis": [20000]*5, "diagnoses.ajcc_pathologic_stage": ["Stage IIA"]*5,
            "demographic.sex_at_birth": ["female"]*5})
        f = pd.DataFrame({"cases.submitter_id": ["p1", "p1"], "follow_ups.days_to_follow_up": [150,250]})
        result = build_os(d, f)
        self.assertEqual(result.loc["p1", "OS_days"], 250)
        self.assertEqual(result.loc["p1", "OS_event"], 0)
        self.assertEqual(result.loc["p1", "OS_exclusion"], "")
        self.assertEqual(result.loc["p2", "OS_exclusion"], "death_without_time")
        self.assertEqual(result.loc["p3", "OS_days"], 400)
        self.assertEqual(result.loc["p4", "OS_exclusion"], "nonpositive_or_missing_time")
        self.assertEqual(len(result), 4)

    def test_stage_parser(self):
        self.assertEqual(stage_number("Stage IIIC"), 3)
        self.assertEqual(stage_number("Stage IV"), 4)
        self.assertTrue(np.isnan(stage_number("Not Reported")))


if __name__ == "__main__":
    unittest.main()

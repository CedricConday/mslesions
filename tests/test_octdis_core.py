from mslesions.octdis.core import dis, optic_nerve


def test_thresholds():
    q = {"qc_pass": True}
    assert optic_nerve({**q, "prnfl_od": 95, "prnfl_os": 89})["optic_nerve"] is True       # 6 um
    assert optic_nerve({**q, "prnfl_od": 95, "prnfl_os": 90})["optic_nerve"] is False      # 5 um: below the consensus 6
    assert optic_nerve({**q, "gcipl_od": 70, "gcipl_os": 66})["optic_nerve"] is True       # 4 um


def test_conditions():
    assert optic_nerve({"prnfl_od": 95, "prnfl_os": 80, "qc_pass": False})["optic_nerve"] is None
    early = optic_nerve({"qc_pass": True, "prnfl_od": 120, "prnfl_os": 95}, last_on_date="2025-01-01", exam_date="2025-02-01")
    assert early["optic_nerve"] is None and "not validated" in early["notes"][0]
    assert optic_nerve({"qc_pass": True, "prnfl_od": 95, "prnfl_os": 80}, other_cause_excluded=False)["optic_nerve"] is None
    assert optic_nerve({}, vep_z=3.0)["optic_nerve"] is True


def test_dis_count():
    d = dis({"optic_nerve": True, "periventricular": True, "juxtacortical_cortical": False, "infratentorial": False})
    assert d["dis_met"] and d["count"] == 2
    d = dis({"periventricular": True, "juxtacortical_cortical": False, "infratentorial": False})
    assert not d["dis_met"] and d["could_change"] and set(d["not_assessed"]) == {"optic_nerve", "spinal_cord"}


def test_review_regressions():
    assert optic_nerve({"prnfl_od": 96, "prnfl_os": 88})["optic_nerve"] is None               # QC not recorded
    assert optic_nerve({"qc_pass": True, "prnfl_od": 90, "prnfl_os": None})["optic_nerve"] is None  # nothing compared
    assert optic_nerve({"qc_pass": True, "bilateral_disease": True, "prnfl_od": 80, "prnfl_os": 81})["optic_nerve"] is None
